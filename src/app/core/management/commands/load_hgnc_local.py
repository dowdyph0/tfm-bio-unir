import pandas as pd
from django.core.management.base import BaseCommand

from core.models import HGNCGene


DEFAULT_PATH = "/app/data/hgnc_complete_set.csv"


class Command(BaseCommand):
    help = "Carga HGNCGene desde fichero local (sin descargar de internet)"

    def add_arguments(self, parser):
        parser.add_argument(
            "--path",
            default=DEFAULT_PATH,
            help="Ruta al fichero TSV/CSV de HGNC",
        )
        parser.add_argument(
            "--batch-size",
            type=int,
            default=2000,
        )

    def handle(self, *args, **options):
        path = options["path"]
        batch_size = options["batch_size"]

        self.stdout.write(f"Cargando HGNC desde: {path}")
        df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        df = df.where(df.notna() & (df != ""), other=None)

        self.stdout.write(f"Filas en fichero: {len(df)}")

        created = 0
        updated = 0

        for start in range(0, len(df), batch_size):
            batch = df.iloc[start : start + batch_size]
            for _, row in batch.iterrows():
                hgnc_id = row.get("hgnc_id")
                if not hgnc_id:
                    continue
                obj, is_new = HGNCGene.objects.update_or_create(
                    hgnc_id=hgnc_id,
                    defaults={
                        "symbol": row.get("symbol") or "",
                        "name": row.get("name") or "",
                        "locus_group": row.get("locus_group") or "",
                        "locus_type": row.get("locus_type") or "",
                        "status": row.get("status") or "",
                        "location": row.get("location") or "",
                        "entrez_id": row.get("entrez_id") or "",
                        "ensembl_gene_id": row.get("ensembl_gene_id") or "",
                        "alias_symbol": row.get("alias_symbol") or "",
                        "prev_symbol": row.get("prev_symbol") or "",
                    },
                )
                if is_new:
                    created += 1
                else:
                    updated += 1

            self.stdout.write(f"  Procesados {min(start + batch_size, len(df))}/{len(df)}...")

        self.stdout.write(self.style.SUCCESS(f"\nHGNC cargado: {created} nuevos, {updated} actualizados."))
