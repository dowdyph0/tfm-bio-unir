import os

import pandas as pd
import requests
from django.core.management.base import BaseCommand

from core.models import TCGACase

XENA_URL = (
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
    "TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix"
)
GDC_CASES_URL = "https://api.gdc.cancer.gov/cases"
BATCH_SIZE = 500


class Command(BaseCommand):
    help = "Pobla TCGACase.molecular_subtype con etiquetas PAM50 de UCSC Xena"

    def add_arguments(self, parser):
        parser.add_argument(
            "--xena-cache",
            default="/app/downloads/brca_clinicalMatrix.tsv",
            help="Ruta local donde guardar/reutilizar el fichero Xena descargado",
        )

    def handle(self, *args, **options):
        xena_cache = options["xena_cache"]

        # 1. recoger uuids de bd
        case_ids = list(TCGACase.objects.values_list("case_id", flat=True))
        self.stdout.write(f"TCGACases en BD: {len(case_ids)}")

        # 2. resolver uuid -> submitter_id via gdc api
        self.stdout.write("Resolviendo barcodes desde GDC...")
        uuid_to_barcode = self._fetch_barcodes(case_ids)
        self.stdout.write(f"  Barcodes obtenidos: {len(uuid_to_barcode)}")

        # 3. descargar / reutilizar clinica xena
        self.stdout.write("Descargando matriz clinica Xena TCGA-BRCA...")
        df_xena = self._load_xena(xena_cache)
        self.stdout.write(f"  Filas Xena: {len(df_xena)}, columnas: {list(df_xena.columns[:6])}")

        if "PAM50Call_RNAseq" not in df_xena.columns:
            self.stdout.write(self.style.ERROR("Columna PAM50Call_RNAseq no encontrada en Xena."))
            self.stdout.write(f"Columnas disponibles: {list(df_xena.columns)}")
            return

        # xena usa barcodes de muestra; se truncan a 12 chars (caso)
        pam50_by_case_barcode: dict[str, str] = {}
        for sample_id, row in df_xena.iterrows():
            case_barcode = str(sample_id)[:12]
            pam50 = row.get("PAM50Call_RNAseq")
            if pd.notna(pam50) and str(pam50).strip():
                pam50_by_case_barcode[case_barcode] = str(pam50).strip()

        self.stdout.write(f"  Muestras con PAM50: {len(pam50_by_case_barcode)}")

        # 4. cruzar y actualizar bd
        updated = 0
        not_found = 0

        cases_to_update = []
        for case in TCGACase.objects.all():
            barcode = uuid_to_barcode.get(case.case_id)
            if not barcode:
                not_found += 1
                continue
            pam50 = pam50_by_case_barcode.get(barcode)
            if pam50:
                case.molecular_subtype = pam50
                cases_to_update.append(case)
                updated += 1
            else:
                not_found += 1

        if cases_to_update:
            TCGACase.objects.bulk_update(cases_to_update, ["molecular_subtype"], batch_size=200)

        self.stdout.write(self.style.SUCCESS("\npopulate_pam50 finalizado"))
        self.stdout.write(f"Actualizados con PAM50 : {updated}")
        self.stdout.write(f"Sin etiqueta PAM50     : {not_found}")

        from collections import Counter
        dist = TCGACase.objects.exclude(molecular_subtype="").values_list(
            "molecular_subtype", flat=True
        )
        for subtype, count in sorted(Counter(dist).items()):
            self.stdout.write(f"  {subtype}: {count}")

    def _fetch_barcodes(self, case_ids: list[str]) -> dict[str, str]:
        # devuelve {uuid: submitter_id} para todos los case_ids
        uuid_to_barcode: dict[str, str] = {}

        for start in range(0, len(case_ids), BATCH_SIZE):
            batch = case_ids[start : start + BATCH_SIZE]
            payload = {
                "filters": {
                    "op": "in",
                    "content": {"field": "case_id", "value": batch},
                },
                "fields": "id,submitter_id",
                "format": "json",
                "size": len(batch),
            }
            r = requests.post(
                GDC_CASES_URL,
                json=payload,
                headers={"Content-Type": "application/json"},
                timeout=30,
            )
            r.raise_for_status()
            hits = r.json()["data"]["hits"]
            for hit in hits:
                uuid_to_barcode[hit["id"]] = hit["submitter_id"]
            self.stdout.write(
                f"  Batch {start}-{start + len(batch)}: {len(hits)} resueltos"
            )

        return uuid_to_barcode

    def _load_xena(self, cache_path: str) -> pd.DataFrame:
        # descarga o reutiliza el fichero clinico de xena
        if not os.path.exists(cache_path):
            self.stdout.write(f"  Descargando desde Xena -> {cache_path}")
            r = requests.get(XENA_URL, timeout=120, stream=True)
            r.raise_for_status()
            with open(cache_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=65536):
                    f.write(chunk)
            self.stdout.write("  Descarga completada.")
        else:
            self.stdout.write(f"  Usando cache: {cache_path}")

        df = pd.read_csv(cache_path, sep="\t", index_col=0, low_memory=False, compression="infer")
        return df
