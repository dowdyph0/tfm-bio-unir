import os

import requests
from django.core.management.base import BaseCommand

from core.models import DownloadStatus, TCGACase, TCGAFile
from core.pipeline import TCGADownloadPipeline

GDC_FILES_URL = "https://api.gdc.cancer.gov/files"
GDC_LINK_BATCH = 400


class Command(BaseCommand):
    help = "Re-importa datos de expresion TCGA desde ficheros locales ya descargados"

    def add_arguments(self, parser):
        parser.add_argument("--project-id", default="TCGA-BRCA")
        parser.add_argument("--downloads-folder", default="/app/downloads")
        parser.add_argument(
            "--num-results",
            type=int,
            default=500,
            help="Numero maximo de ficheros a enqueue via GDC API",
        )
        parser.add_argument(
            "--data-types",
            nargs="+",
            default=["Gene Expression Quantification"],
        )

    def handle(self, *args, **options):
        downloads_folder = options["downloads_folder"]
        data_types = options["data_types"]

        pipeline = TCGADownloadPipeline()

        # 1. indice de ficheros locales (file_id -> ruta)
        self.stdout.write("Indexando ficheros locales en disco...")
        local_index: dict[str, str] = {}
        if os.path.isdir(downloads_folder):
            for fname in os.listdir(downloads_folder):
                # los ficheros se llaman {file_id}.rna_seq.augmented_star_gene_counts.tsv
                file_id_candidate = fname.split(".")[0]
                if file_id_candidate:
                    local_index[file_id_candidate] = os.path.join(downloads_folder, fname)
        self.stdout.write(f"  {len(local_index)} ficheros encontrados en disco.")

        # 2. crear TCGAFile desde los ficheros en disco (sin api)
        EXPR_SUFFIX = ".rna_seq.augmented_star_gene_counts.tsv"
        created_from_disk = 0
        for fname, fpath in local_index.items():
            filename = os.path.basename(fpath)
            if not filename.endswith(EXPR_SUFFIX):
                continue
            _, created = TCGAFile.objects.get_or_create(
                file_id=fname,
                defaults={
                    "file_name": filename,
                    "data_type": "Gene Expression Quantification",
                    "experimental_strategy": "RNA-Seq",
                    "file_size": os.path.getsize(fpath),
                    "local_path": fpath,
                },
            )
            if created:
                created_from_disk += 1

        self.stdout.write(f"  TCGAFile creados desde disco: {created_from_disk}")

        # 3. parsear todos los ficheros de expresion con ruta local
        qs = TCGAFile.objects.filter(data_type__in=data_types)
        total = qs.count()
        self.stdout.write(f"\nTotal TCGAFile en BD: {total}")

        found = 0
        parsed = 0
        failed = 0
        skipped_no_file = 0

        for tcga_file in qs.iterator():
            local_path = local_index.get(tcga_file.file_id)
            if not local_path:
                skipped_no_file += 1
                self.stdout.write(
                    self.style.WARNING(f"  [SKIP] {tcga_file.file_id} - no encontrado en disco")
                )
                continue

            found += 1
            tcga_file.local_path = local_path
            tcga_file.save(update_fields=["local_path"])

            status_obj, _ = DownloadStatus.objects.get_or_create(tcga_file=tcga_file)
            status_obj.local_path = local_path
            status_obj.downloaded_bytes = os.path.getsize(local_path)
            status_obj.status = DownloadStatus.STATUS_DOWNLOADED
            status_obj.save(update_fields=["local_path", "downloaded_bytes", "status", "updated_at"])

            ok = pipeline._parse_file(tcga_file, status_obj)
            if ok:
                parsed += 1
                self.stdout.write(f"  [OK]   {tcga_file.file_name}")
            else:
                failed += 1
                self.stdout.write(
                    self.style.ERROR(f"  [FAIL] {tcga_file.file_name}: {status_obj.error_message[:120]}")
                )

        self.stdout.write(self.style.SUCCESS("\nReimport finalizado"))
        self.stdout.write(f"TCGAFile creados desde disco: {created_from_disk}")
        self.stdout.write(f"Ficheros en BD     : {total}")
        self.stdout.write(f"Encontrados en disco: {found}")
        self.stdout.write(f"Parseados OK       : {parsed}")
        self.stdout.write(f"Fallidos           : {failed}")
        self.stdout.write(f"Sin fichero local  : {skipped_no_file}")

        # 4. vincular TCGAFile <-> TCGACase via gdc api
        self.stdout.write("\nVinculando TCGAFile <-> TCGACase via GDC API...")
        ThroughModel = TCGAFile.cases.through
        case_by_uuid = {c.case_id: c for c in TCGACase.objects.all()}
        # solo ficheros con local_path
        all_file_ids = list(
            TCGAFile.objects
            .filter(data_type__in=data_types, local_path__isnull=False)
            .exclude(local_path="")
            .values_list("file_id", flat=True)
        )
        m2m_linked = 0
        for i in range(0, len(all_file_ids), GDC_LINK_BATCH):
            batch = all_file_ids[i:i + GDC_LINK_BATCH]
            payload = {
                "filters": {"op": "in", "content": {"field": "file_id", "value": batch}},
                "fields": "file_id,cases.case_id",
                "size": GDC_LINK_BATCH,
                "format": "json",
            }
            try:
                resp = requests.post(GDC_FILES_URL, json=payload, timeout=60)
                resp.raise_for_status()
            except Exception as exc:
                self.stdout.write(self.style.WARNING(f"  [WARN] Lote {i//GDC_LINK_BATCH+1} fallo: {exc}"))
                continue
            hits = resp.json()["data"]["hits"]
            file_by_uuid = {
                f.file_id: f for f in TCGAFile.objects.filter(file_id__in=batch)
            }
            to_create = []
            for hit in hits:
                tcga_file = file_by_uuid.get(hit["file_id"])
                if not tcga_file:
                    continue
                for case_info in hit.get("cases", []):
                    case_obj = case_by_uuid.get(case_info.get("case_id", ""))
                    if case_obj:
                        to_create.append(ThroughModel(
                            tcgafile_id=tcga_file.pk,
                            tcgacase_id=case_obj.pk,
                        ))
            if to_create:
                ThroughModel.objects.bulk_create(to_create, ignore_conflicts=True)
                m2m_linked += len(to_create)
            self.stdout.write(f"  Lote {i//GDC_LINK_BATCH+1}: {len(hits)} ficheros -> {len(to_create)} relaciones")

        self.stdout.write(self.style.SUCCESS(f"M2M creadas: {m2m_linked} | Total en BD: {ThroughModel.objects.count()}"))
