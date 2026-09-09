import logging

from django.core.management.base import BaseCommand

logger = logging.getLogger(__name__)

PAM50_LABELS = {"Luminal A", "Luminal B", "HER2-enriched", "Basal-like", "Normal-like"}

# mapeo de etiquetas geo a nombres cortos pam50
PAM50_NORMALIZE = {
    "luminal a": "LumA",
    "luma": "LumA",
    "luminal b": "LumB",
    "lumb": "LumB",
    "her2-enriched": "Her2",
    "her2": "Her2",
    "basal-like": "Basal",
    "basal": "Basal",
    "normal-like": "Normal",
    "normal": "Normal",
}


class Command(BaseCommand):
    help = "Descarga una serie GEO, mapea genes via HGNC y calcula z-scores"

    def add_arguments(self, parser):
        parser.add_argument("--gse", default="GSE25055",
                            help="Accession de la serie GEO (default: GSE25055)")
        parser.add_argument("--destdir", default="/tmp/geo",
                            help="Directorio de cache para ficheros GEO descargados")
        parser.add_argument("--limit-samples", type=int, default=0,
                            help="Limitar numero de muestras (0 = todas, util para pruebas)")

    def handle(self, *args, **options):
        import os
        import GEOparse
        from core.models import GEOSeries, GEOSample, GEOExpressionRecord, HGNCGene

        gse_id = options["gse"]
        destdir = options["destdir"]
        limit = options["limit_samples"]
        os.makedirs(destdir, exist_ok=True)

        self.stdout.write(f"[sync_geo] Descargando {gse_id} en {destdir} ...")
        gse = GEOparse.get_GEO(geo=gse_id, destdir=destdir, silent=False)

        # 1. geoseries en bd
        series_obj, created = GEOSeries.objects.get_or_create(
            accession=gse_id,
            defaults={
                "title": gse.metadata.get("title", [""])[0],
                "gds_type": gse.metadata.get("type", [""])[0],
                "taxon": gse.metadata.get("organism", [""])[0],
                "n_samples": len(gse.gsms),
                "platform_id": gse.metadata.get("platform_id", [""])[0],
            }
        )
        self.stdout.write(f"  GEOSeries {'creada' if created else 'ya existe'}: {gse_id}")

        # 2. plataforma: probe_id -> entrez_id
        platform_id = gse.metadata.get("platform_id", [""])[0]
        self.stdout.write(f"  Descargando plataforma {platform_id} ...")
        gpl = GEOparse.get_GEO(geo=platform_id, destdir=destdir, silent=True)

        entrez_col = None
        for col in gpl.table.columns:
            if "entrez" in col.lower() or "gene_id" in col.lower() or col == "ENTREZ_GENE_ID":
                entrez_col = col
                break

        if entrez_col is None:
            # fallback: columna con numeros que parezcan entrez ids
            for col in gpl.table.columns:
                if col.upper() not in ("ID", "GB_ACC", "SPOT_ID", "SEQUENCE"):
                    entrez_col = col
                    break

        self.stdout.write(f"  Columna Entrez en GPL: {entrez_col}")

        # probe_id -> entrez_id (puede haber varios separados por ///)
        probe_to_entrez = {}
        for _, row in gpl.table.iterrows():
            probe = str(row.get("ID", row.name))
            entrez_raw = str(row.get(entrez_col, "")).strip()
            if not entrez_raw or entrez_raw in ("nan", "---", ""):
                continue
            # solo el primero si hay varios
            entrez = entrez_raw.split("///")[0].strip()
            probe_to_entrez[probe] = entrez

        self.stdout.write(f"  Probes con Entrez ID: {len(probe_to_entrez):,}")

        # 3. hgnc: entrez_id -> hgnc gene
        entrez_ids = set(probe_to_entrez.values())
        hgnc_by_entrez = {
            g.entrez_id: g
            for g in HGNCGene.objects.filter(entrez_id__in=entrez_ids)
        }
        self.stdout.write(f"  Entrez IDs con mapeo HGNC: {len(hgnc_by_entrez):,}")

        # 4. muestras: gsm -> geosample + extraer pam50
        gsm_ids = list(gse.gsms.keys())
        if limit:
            gsm_ids = gsm_ids[:limit]
            self.stdout.write(f"  [--limit-samples] Usando {limit} muestras")

        self.stdout.write(f"  Procesando {len(gsm_ids)} muestras ...")

        sample_objs = {}
        pam50_by_gsm = {}

        for gsm_id in gsm_ids:
            gsm = gse.gsms[gsm_id]
            chars = gsm.metadata.get("characteristics_ch1", [])
            if isinstance(chars, str):
                chars = [chars]

            clinical = {}
            pam50 = ""
            for c in chars:
                if ":" in c:
                    k, v = c.split(":", 1)
                    k_low = k.strip().lower()
                    v_str = v.strip()
                    clinical[k_low] = v_str
                    # detectar campo pam50 (distintos nombres segun serie)
                    if k_low in ("pam50_class", "pam50", "subtype", "molecular subtype",
                                 "molecular_subtype", "pam50 subtype"):
                        norm = PAM50_NORMALIZE.get(v_str.lower(), "")
                        if norm:
                            pam50 = norm

            if pam50:
                clinical["pam50"] = pam50
                pam50_by_gsm[gsm_id] = pam50

            obj, _ = GEOSample.objects.get_or_create(
                gsm_id=gsm_id,
                defaults={
                    "series": series_obj,
                    "title": gsm.metadata.get("title", [""])[0] if gsm.metadata.get("title") else "",
                    "source_name": gsm.metadata.get("source_name_ch1", [""])[0] if gsm.metadata.get("source_name_ch1") else "",
                    "organism": gsm.metadata.get("organism_ch1", [""])[0] if gsm.metadata.get("organism_ch1") else "",
                    "platform_id": gsm.metadata.get("platform_id", [""])[0] if gsm.metadata.get("platform_id") else "",
                    "clinical_metadata": clinical,
                }
            )
            # actualizar pam50 si ya existia el objeto
            if pam50 and obj.clinical_metadata.get("pam50") != pam50:
                obj.clinical_metadata["pam50"] = pam50
                obj.save(update_fields=["clinical_metadata"])

            sample_objs[gsm_id] = obj

        self.stdout.write(f"  Muestras con etiqueta PAM50: {len(pam50_by_gsm)}/{len(gsm_ids)}")

        # 5. matriz de expresion: pivot muestras x probes
        self.stdout.write("  Construyendo matriz de expresion ...")
        import pandas as pd

        frames = []
        for gsm_id in gsm_ids:
            gsm = gse.gsms[gsm_id]
            if gsm.table is None or gsm.table.empty:
                continue
            # columna de valor: VALUE o la segunda columna
            val_col = "VALUE" if "VALUE" in gsm.table.columns else gsm.table.columns[1]
            s = gsm.table.set_index("ID_REF")[val_col].rename(gsm_id)
            frames.append(s)

        if not frames:
            self.stderr.write("ERROR: no se pudo construir la matriz (tablas vacias)")
            return

        df_expr = pd.concat(frames, axis=1)  # probes x muestras
        self.stdout.write(f"  Matriz bruta: {df_expr.shape[0]:,} probes x {df_expr.shape[1]} muestras")

        # convertir a numerico
        df_expr = df_expr.apply(pd.to_numeric, errors="coerce")

        # 6. probe_id -> gene_symbol (via entrez -> hgnc)
        probe_to_symbol = {}
        probe_to_hgnc = {}
        for probe in df_expr.index:
            entrez = probe_to_entrez.get(str(probe))
            if not entrez:
                continue
            hgnc = hgnc_by_entrez.get(entrez)
            if hgnc:
                probe_to_symbol[probe] = hgnc.symbol
                probe_to_hgnc[probe] = hgnc

        df_expr = df_expr.loc[df_expr.index.isin(probe_to_symbol)]
        df_expr.index = df_expr.index.map(probe_to_symbol)
        # varios probes por gen -> media
        df_expr = df_expr.groupby(level=0).mean()
        self.stdout.write(f"  Genes unicos mapeados a HGNC: {df_expr.shape[0]:,}")

        # 7. z-score por gen (sobre todas las muestras de la serie)
        gene_mean = df_expr.mean(axis=1)
        gene_std = df_expr.std(axis=1).replace(0, 1)
        df_zscore = df_expr.subtract(gene_mean, axis=0).divide(gene_std, axis=0)

        # 8. bulk create GEOExpressionRecord
        # borrar registros previos de esta serie para permitir re-sync
        old_ids = list(GEOSample.objects.filter(series=series_obj).values_list("id", flat=True))
        if old_ids:
            deleted, _ = GEOExpressionRecord.objects.filter(geo_sample_id__in=old_ids).delete()
            if deleted:
                self.stdout.write(f"  Borrados {deleted} registros previos")

        BATCH = 5000
        total = 0
        records = []

        symbol_to_hgnc = {
            probe_to_symbol[p]: probe_to_hgnc[p]
            for p in probe_to_hgnc
        }

        for gene_symbol in df_expr.index:
            hgnc = symbol_to_hgnc.get(gene_symbol)
            for gsm_id in df_expr.columns:
                sample_obj = sample_objs.get(gsm_id)
                if sample_obj is None:
                    continue
                val = df_expr.at[gene_symbol, gsm_id]
                zsc = df_zscore.at[gene_symbol, gsm_id]
                if pd.isna(val):
                    continue
                records.append(GEOExpressionRecord(
                    geo_sample=sample_obj,
                    hgnc_gene=hgnc,
                    gene_symbol=gene_symbol,
                    value=float(val),
                    zscore=float(zsc) if not pd.isna(zsc) else None,
                ))
                if len(records) >= BATCH:
                    GEOExpressionRecord.objects.bulk_create(records, ignore_conflicts=True)
                    total += len(records)
                    records = []
                    self.stdout.write(f"  ... {total:,} registros insertados")

        if records:
            GEOExpressionRecord.objects.bulk_create(records, ignore_conflicts=True)
            total += len(records)

        self.stdout.write(self.style.SUCCESS(
            f"\n[sync_geo] Completado: {total:,} GEOExpressionRecord insertados "
            f"({df_expr.shape[0]:,} genes x {len(gsm_ids)} muestras, "
            f"{len(pam50_by_gsm)} con etiqueta PAM50)"
        ))
