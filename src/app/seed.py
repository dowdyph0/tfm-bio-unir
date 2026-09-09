import os
import sys

sys.path.insert(0, "/app")

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'project.settings')
import django
django.setup()

from core.connectors import TCGAConnector, GEOConnector, HGNCConnector
from core.models import TCGAProject, TCGACase, TCGAFile, GEOSample, HGNCGene, GEOSeries


def seed_tcga():
    print("[TCGA] Descargando proyectos...")
    tcga = TCGAConnector()

    # list_projects devuelve project_id, name, disease_type, primary_site,
    # summary.case_count y summary.file_count
    df_proj = tcga.list_projects(keyword="breast")

    for _, row in df_proj.iterrows():
        obj, created = TCGAProject.objects.get_or_create(
            project_id=row["project_id"],
            defaults={
                "name": row.get("name", ""),
                "primary_site": row.get("primary_site", []),
                "disease_type": row.get("disease_type", []),
                "case_count": row.get("summary.case_count", 0),
                "file_count": row.get("summary.file_count", 0),
            }
        )
        print(f"  {'Creado' if created else 'Ya existe'}: {obj.project_id}")

    first_project_id = df_proj["project_id"].iloc[0]
    project_obj = TCGAProject.objects.get(project_id=first_project_id)

    print(f"[TCGA] Descargando casos del proyecto {first_project_id}...")
    df_cases = tcga.get_clinical_summary(project_id=first_project_id, num_results=3)

    for _, row in df_cases.iterrows():
        case_obj, created = TCGACase.objects.get_or_create(
            case_id=row["case_id"],
            defaults={
                "project": project_obj,
                "primary_site": row.get("primary_site", ""),
                "disease_type": row.get("disease_type", ""),
                "gender": row.get("demographic.gender", ""),
                "age_at_index": row.get("demographic.age_at_index") or None,
            }
        )
        print(f"  {'Creado' if created else 'Ya existe'}: {case_obj.case_id}")

        # list_files filtra por project_id; no hay metodo por case_id,
        # asi que se asocian los primeros ficheros del proyecto al caso
        df_files = tcga.list_files(project_id=first_project_id, num_results=2)
        for _, frow in df_files.iterrows():
            file_obj, _ = TCGAFile.objects.get_or_create(
                file_id=frow["file_id"],
                defaults={
                    "file_name": frow.get("file_name", ""),
                    "data_type": frow.get("data_type", ""),
                    "experimental_strategy": frow.get("experimental_strategy", ""),
                    "file_size": frow.get("file_size", 0),
                }
            )
            file_obj.cases.add(case_obj)

    print("[TCGA] Seed completado.\n")


def seed_geo():
    print("[GEO] Descargando metadatos de serie de prueba...")
    geo = GEOConnector()

    studies = geo.search(
        "breast cancer[title] AND expression profiling by array[DataSet Type]",
        num_results=1
    )
    if not studies:
        print("  No se encontraron series GEO.")
        return

    study = studies[0]
    gse_id = study["accession"]
    print(f"  Serie seleccionada: {gse_id} - {study['title']}")

    series_obj, created = GEOSeries.objects.get_or_create(
        accession=gse_id,
        defaults={
            "title": study.get("title", ""),
            "gds_type": study.get("gdstype", ""),
            "taxon": study.get("taxon", ""),
            "n_samples": study.get("nsamples", 0),
            "publication_date": study.get("pdat", ""),
        }
    )
    print(f"  GEOSeries {'creada' if created else 'ya existe'}: {series_obj.accession}")

    metadata = geo.get_series_metadata(gse_id)
    df_samples = metadata["samples_df"].head(3)

    def first(val):
        if isinstance(val, list):
            return val[0] if val else ""
        return val or ""

    for _, row in df_samples.iterrows():
        gsm_id = row.get("geo_accession") or row.name

        obj, created = GEOSample.objects.get_or_create(
            gsm_id=gsm_id,
            defaults={
                "series": series_obj,
                "title": first(row.get("title")),
                "source_name": first(row.get("source_name_ch1")),
                "organism": first(row.get("organism_ch1")),
                "platform_id": first(row.get("platform_id")),
                "clinical_metadata": row.get("characteristics_ch1") or {},
            }
        )
        print(f"  GEOSample {'creada' if created else 'ya existe'}: {gsm_id}")

    print("[GEO] Seed completado.\n")


def seed_hgnc():
    print("[HGNC] Descargando y cargando genes...")
    hgnc = HGNCConnector()
    df = hgnc.download_complete_set(dest_path="/tmp/hgnc_complete_set.txt")

    count = 0
    for _, row in df.iterrows():
        _, created = HGNCGene.objects.get_or_create(
            hgnc_id=row["hgnc_id"],
            defaults={
                "symbol": row.get("symbol"),
                "name": row.get("name"),
                "entrez_id": row.get("entrez_id"),
                "ensembl_gene_id": row.get("ensembl_gene_id"),
            }
        )
        if created:
            count += 1
        if count >= 50:  # limite seed de prueba
            break
    print(f"[HGNC] {count} genes cargados.\n")


if __name__ == "__main__":
    seed_hgnc()
    seed_tcga()
    seed_geo()