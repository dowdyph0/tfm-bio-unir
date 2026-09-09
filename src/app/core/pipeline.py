from __future__ import annotations

import os
import math

from django.db import transaction
from django.utils import timezone

from .connectors import TCGAConnector
from .models import (
    DownloadStatus,
    HGNCGene,
    TCGAExpressionRecord,
    TCGAFile,
    TCGAMutation,
)
from .parsers import parse_tcga_file


class TCGADownloadPipeline:
    # descarga y parsea ficheros tcga con seguimiento en DownloadStatus

    def __init__(self, connector: TCGAConnector | None = None):
        self.connector = connector or TCGAConnector()

    def run_for_existing_files(
        self,
        queryset,
        dest_folder: str,
        parse_downloaded: bool = True,
    ) -> dict:
        os.makedirs(dest_folder, exist_ok=True)

        summary = {
            "total": 0,
            "downloaded": 0,
            "parsed": 0,
            "failed": 0,
            "skipped": 0,
        }

        for tcga_file in queryset:
            summary["total"] += 1
            status_obj, _ = DownloadStatus.objects.get_or_create(tcga_file=tcga_file)

            if tcga_file.local_path and os.path.exists(tcga_file.local_path):
                status_obj.local_path = tcga_file.local_path
                status_obj.downloaded_bytes = os.path.getsize(tcga_file.local_path)
                if status_obj.status not in [DownloadStatus.STATUS_PARSED, DownloadStatus.STATUS_FAILED]:
                    status_obj.status = DownloadStatus.STATUS_DOWNLOADED
                    status_obj.save(update_fields=["local_path", "downloaded_bytes", "status", "updated_at"])
                summary["skipped"] += 1
                if parse_downloaded:
                    parsed_ok = self._parse_file(tcga_file, status_obj)
                    if parsed_ok:
                        summary["parsed"] += 1
                    else:
                        summary["failed"] += 1
                continue

            ok = self._download_and_optionally_parse(tcga_file, status_obj, dest_folder, parse_downloaded)
            if ok == "parsed":
                summary["downloaded"] += 1
                summary["parsed"] += 1
            elif ok == "downloaded":
                summary["downloaded"] += 1
            else:
                summary["failed"] += 1

        return summary

    def _download_and_optionally_parse(
        self,
        tcga_file: TCGAFile,
        status_obj: DownloadStatus,
        dest_folder: str,
        parse_downloaded: bool,
    ) -> str:
        status_obj.status = DownloadStatus.STATUS_DOWNLOADING
        status_obj.started_at = timezone.now()
        status_obj.error_message = ""
        status_obj.save(update_fields=["status", "started_at", "error_message", "updated_at"])

        try:
            result = self.connector.download_file(file_id=tcga_file.file_id, dest_folder=dest_folder)
            tcga_file.local_path = result["local_path"]
            tcga_file.save(update_fields=["local_path"])

            status_obj.local_path = result["local_path"]
            status_obj.downloaded_bytes = int(result.get("downloaded_bytes", 0))
            status_obj.finished_at = timezone.now()
            status_obj.status = DownloadStatus.STATUS_DOWNLOADED
            status_obj.save(
                update_fields=[
                    "local_path",
                    "downloaded_bytes",
                    "finished_at",
                    "status",
                    "updated_at",
                ]
            )

            if parse_downloaded:
                parsed_ok = self._parse_file(tcga_file, status_obj)
                if parsed_ok:
                    return "parsed"
                return "failed"

            return "downloaded"
        except Exception as exc:
            status_obj.status = DownloadStatus.STATUS_FAILED
            status_obj.error_message = str(exc)
            status_obj.finished_at = timezone.now()
            status_obj.save(update_fields=["status", "error_message", "finished_at", "updated_at"])
            return "failed"

    def _parse_file(self, tcga_file: TCGAFile, status_obj: DownloadStatus) -> bool:
        try:
            parser_used, parse_summary, df = parse_tcga_file(
                file_path=tcga_file.local_path,
                file_name=tcga_file.file_name,
                data_type=tcga_file.data_type,
                include_dataframe=True,
            )
            persisted_rows = self._persist_parsed_data(tcga_file=tcga_file, parser_used=parser_used, df=df)
            parse_summary["persisted_rows"] = persisted_rows

            status_obj.parser_used = parser_used
            status_obj.parse_summary = parse_summary
            status_obj.status = DownloadStatus.STATUS_PARSED if parser_used != "none" else DownloadStatus.STATUS_DOWNLOADED
            status_obj.save(update_fields=["parser_used", "parse_summary", "status", "updated_at"])
            return True
        except Exception as exc:
            status_obj.status = DownloadStatus.STATUS_FAILED
            status_obj.error_message = f"Parser error: {exc}"
            status_obj.finished_at = timezone.now()
            status_obj.save(update_fields=["status", "error_message", "finished_at", "updated_at"])
            return False

    def _persist_parsed_data(self, tcga_file: TCGAFile, parser_used: str, df) -> int:
        if df is None:
            return 0
        if parser_used == "maf":
            return self._persist_mutations(tcga_file, df)
        if parser_used == "expression":
            return self._persist_expression(tcga_file, df)
        return 0

    @transaction.atomic
    def _persist_mutations(self, tcga_file: TCGAFile, df) -> int:
        TCGAMutation.objects.filter(tcga_file=tcga_file).delete()

        hgnc_by_symbol = {g.symbol: g for g in HGNCGene.objects.filter(symbol__in=df.get("Hugo_Symbol", []).dropna().unique())} if "Hugo_Symbol" in df.columns else {}

        objects = []
        for _, row in df.iterrows():
            symbol = str(row.get("Hugo_Symbol") or "").strip()
            start_position = self._safe_int(row.get("Start_Position"))
            end_position = self._safe_int(row.get("End_Position"))

            objects.append(
                TCGAMutation(
                    tcga_file=tcga_file,
                    hgnc_gene=hgnc_by_symbol.get(symbol),
                    hugo_symbol=symbol,
                    entrez_gene_id=self._safe_str(row.get("Entrez_Gene_Id")),
                    variant_classification=self._safe_str(row.get("Variant_Classification")),
                    variant_type=self._safe_str(row.get("Variant_Type")),
                    chromosome=self._safe_str(row.get("Chromosome")),
                    start_position=start_position,
                    end_position=end_position,
                    reference_allele=self._safe_str(row.get("Reference_Allele")),
                    tumor_seq_allele2=self._safe_str(row.get("Tumor_Seq_Allele2")),
                    tumor_sample_barcode=self._safe_str(row.get("Tumor_Sample_Barcode")),
                    raw_row={k: self._json_safe(v) for k, v in row.items()},
                )
            )

        TCGAMutation.objects.bulk_create(objects, batch_size=1000)
        return len(objects)

    @transaction.atomic
    def _persist_expression(self, tcga_file: TCGAFile, df) -> int:
        TCGAExpressionRecord.objects.filter(tcga_file=tcga_file).delete()

        ensembl_ids = []
        if "gene_id" in df.columns:
            for raw_gene_id in df["gene_id"].dropna().astype(str).tolist():
                ensembl_ids.append(raw_gene_id.split(".")[0])
        hgnc_by_ensembl = {
            g.ensembl_gene_id: g
            for g in HGNCGene.objects.filter(ensembl_gene_id__in=ensembl_ids)
            if g.ensembl_gene_id
        }

        objects = []
        for _, row in df.iterrows():
            gene_id = self._safe_str(row.get("gene_id"))
            gene_id_versionless = gene_id.split(".")[0] if gene_id else ""

            objects.append(
                TCGAExpressionRecord(
                    tcga_file=tcga_file,
                    hgnc_gene=hgnc_by_ensembl.get(gene_id_versionless),
                    gene_id=gene_id,
                    gene_id_versionless=gene_id_versionless,
                    gene_name=self._safe_str(row.get("gene_name")),
                    gene_type=self._safe_str(row.get("gene_type")),
                    unstranded=self._safe_float(row.get("unstranded")),
                    stranded_first=self._safe_float(row.get("stranded_first")),
                    stranded_second=self._safe_float(row.get("stranded_second")),
                    tpm_unstranded=self._safe_float(row.get("tpm_unstranded")),
                    fpkm_unstranded=self._safe_float(row.get("fpkm_unstranded")),
                    fpkm_uq_unstranded=self._safe_float(row.get("fpkm_uq_unstranded")),
                    raw_row={k: self._json_safe(v) for k, v in row.items()},
                )
            )

        TCGAExpressionRecord.objects.bulk_create(objects, batch_size=1000)
        return len(objects)

    @staticmethod
    def _safe_int(value):
        try:
            if value is None or value == "":
                return None
            return int(value)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _safe_float(value):
        try:
            if value is None or value == "":
                return None
            normalized = float(value)
            if math.isnan(normalized):
                return None
            return normalized
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _safe_str(value):
        if value is None:
            return ""
        if isinstance(value, float) and math.isnan(value):
            return ""
        return str(value).strip()

    @staticmethod
    def _json_safe(value):
        if value is None:
            return None
        if isinstance(value, dict):
            return {k: TCGADownloadPipeline._json_safe(v) for k, v in value.items()}
        if isinstance(value, list):
            return [TCGADownloadPipeline._json_safe(v) for v in value]
        if isinstance(value, float) and math.isnan(value):
            return None
        if hasattr(value, "item"):
            try:
                normalized = value.item()
                if isinstance(normalized, float) and math.isnan(normalized):
                    return None
                return normalized
            except Exception:
                return str(value)
        return value

    def enqueue_project_files(
        self,
        project_id: str,
        data_type: str,
        num_results: int = 10,
    ) -> int:
        # importa metadatos de tcga files a bd para luego descargarlos en pipeline
        df_files = self.connector.list_files(project_id=project_id, data_type=data_type, num_results=num_results)
        created_or_updated = 0

        for _, row in df_files.iterrows():
            tcga_file, _ = TCGAFile.objects.update_or_create(
                file_id=row["file_id"],
                defaults={
                    "file_name": row.get("file_name", ""),
                    "data_type": row.get("data_type", ""),
                    "experimental_strategy": row.get("experimental_strategy", ""),
                    "file_size": row.get("file_size", 0) or 0,
                },
            )
            DownloadStatus.objects.get_or_create(tcga_file=tcga_file)
            created_or_updated += 1

        return created_or_updated
