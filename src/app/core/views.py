import csv
import io
import logging
import os
import zipfile
from datetime import datetime, timezone

import pandas as pd
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.db.models import Q
from django.http import HttpResponse
from django.shortcuts import render
from django.views import View

from .models import HGNCGene, TCGAExpressionRecord, TCGAMutation

logger = logging.getLogger(__name__)


class GeneIntegratedView(LoginRequiredMixin, UserPassesTestMixin, View):
	login_url = "/admin/login/"
	template_name = "core/gene_integrated_view.html"

	def test_func(self):
		return self.request.user.is_staff

	MUTATION_EXPORT_FIELDS = [
		"hugo_symbol",
		"entrez_gene_id",
		"variant_classification",
		"variant_type",
		"chromosome",
		"start_position",
		"end_position",
		"tumor_sample_barcode",
		"tcga_file__file_id",
		"tcga_file__file_name",
	]

	EXPRESSION_EXPORT_FIELDS = [
		"gene_id",
		"gene_id_versionless",
		"gene_name",
		"gene_type",
		"tpm_unstranded",
		"fpkm_unstranded",
		"fpkm_uq_unstranded",
		"tcga_file__file_id",
		"tcga_file__file_name",
	]

	def get(self, request):
		symbol = request.GET.get("symbol", "").strip()
		selected_project = request.GET.get("project_id", "").strip()
		selected_variant = request.GET.get("variant_classification", "").strip()
		min_tpm = self._parse_float(request.GET.get("min_tpm", ""))
		min_fpkm_uq = self._parse_float(request.GET.get("min_fpkm_uq", ""))
		row_limit = request.GET.get("limit", "100").strip()
		export_kind = request.GET.get("export", "").strip()

		logger.info(
			"Gene integrated view request user=%s method=%s symbol=%s project=%s variant=%s min_tpm=%s min_fpkm_uq=%s limit=%s export=%s",
			request.user.get_username() if request.user.is_authenticated else "anonymous",
			request.method,
			symbol,
			selected_project,
			selected_variant,
			min_tpm,
			min_fpkm_uq,
			row_limit,
			export_kind,
		)

		gene = None
		mutations = []
		expressions = []
		available_raw_count = 0
		missing_raw_paths = []
		project_options = []
		variant_options = []
		total_mutations = 0
		total_expressions = 0
		exported_at = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

		if symbol:
			gene = HGNCGene.objects.filter(symbol__iexact=symbol).first()
			logger.info("Gene lookup result symbol=%s found_hgnc=%s", symbol, bool(gene))

			mutation_filter = Q(hugo_symbol__iexact=symbol)
			expression_filter = Q(gene_name__iexact=symbol)

			if gene:
				mutation_filter = mutation_filter | Q(hgnc_gene=gene)
				expression_filter = expression_filter | Q(hgnc_gene=gene)
				if gene.ensembl_gene_id:
					expression_filter = expression_filter | Q(gene_id_versionless=gene.ensembl_gene_id)

			mutation_base_qs = (
				TCGAMutation.objects.filter(mutation_filter)
				.select_related("tcga_file")
				.order_by("tcga_file__file_name")
				.distinct()
			)
			expression_base_qs = (
				TCGAExpressionRecord.objects.filter(expression_filter)
				.select_related("tcga_file")
				.order_by("tcga_file__file_name")
				.distinct()
			)

			raw_project_values = list(
				mutation_base_qs.values_list("tcga_file__cases__project__project_id", flat=True)
			) + list(
				expression_base_qs.values_list("tcga_file__cases__project__project_id", flat=True)
			)
			project_options = sorted(
				{
					str(p).strip()
					for p in raw_project_values
					if p is not None and str(p).strip()
				}
			)

			variant_options = list(
				mutation_base_qs.exclude(variant_classification__isnull=True)
				.exclude(variant_classification="")
				.values_list("variant_classification", flat=True)
				.distinct()
				.order_by("variant_classification")
			)

			filtered_mutations = mutation_base_qs
			filtered_expressions = expression_base_qs

			if selected_project:
				filtered_mutations = filtered_mutations.filter(
					tcga_file__cases__project__project_id=selected_project
				).distinct()
				filtered_expressions = filtered_expressions.filter(
					tcga_file__cases__project__project_id=selected_project
				).distinct()

			if selected_variant:
				filtered_mutations = filtered_mutations.filter(variant_classification=selected_variant)

			if min_tpm is not None:
				filtered_expressions = filtered_expressions.filter(tpm_unstranded__gte=min_tpm)

			if min_fpkm_uq is not None:
				filtered_expressions = filtered_expressions.filter(fpkm_uq_unstranded__gte=min_fpkm_uq)

			total_mutations = filtered_mutations.count()
			total_expressions = filtered_expressions.count()

			filtered_mutations, filtered_expressions = self._apply_limit(
				filtered_mutations,
				filtered_expressions,
				row_limit,
			)

			mutations = list(filtered_mutations)
			expressions = list(filtered_expressions)
			self._enrich_rows_with_links(mutations, expressions)

			raw_entries, missing_raw_paths = self._build_raw_entries(
				filtered_mutations,
				filtered_expressions,
			)
			available_raw_count = len(raw_entries)

			if export_kind:
				response = self._build_export_response(
					export_kind=export_kind,
					gene=gene,
					symbol=symbol,
					exported_at=exported_at,
					mutations_qs=filtered_mutations,
					expressions_qs=filtered_expressions,
					raw_entries=raw_entries,
					missing_raw_paths=missing_raw_paths,
					request=request,
				)
				if response is not None:
					return response

			logger.info(
				"Gene data loaded symbol=%s mutations=%d expressions=%d",
				symbol,
				total_mutations,
				total_expressions,
			)

		context = dict(
			symbol=symbol,
			gene=gene,
			gene_links=self._build_gene_links(gene=gene, symbol=symbol, mutations=mutations, expressions=expressions),
			mutations=mutations,
			expressions=expressions,
			selected_project=selected_project,
			selected_variant=selected_variant,
			min_tpm=request.GET.get("min_tpm", "").strip(),
			min_fpkm_uq=request.GET.get("min_fpkm_uq", "").strip(),
			row_limit=row_limit,
			project_options=project_options,
			variant_options=variant_options,
			total_mutations=total_mutations,
			total_expressions=total_expressions,
			available_raw_count=available_raw_count,
			missing_raw_paths=missing_raw_paths,
		)
		return render(request, self.template_name, context)

	def _normalize_ensembl_gene_id(self, value):
		text = (value or "").strip()
		if not text:
			return ""
		return text.split(".")[0]

	def _normalize_entrez_id(self, value):
		text = (value or "").strip()
		return text if text.isdigit() else ""

	def _extract_ensembl_from_expression(self, expression_row):
		if not expression_row:
			return ""
		if getattr(expression_row, "gene_id_versionless", ""):
			return self._normalize_ensembl_gene_id(expression_row.gene_id_versionless)
		return self._normalize_ensembl_gene_id(getattr(expression_row, "gene_id", ""))

	def _build_gene_links(self, gene, symbol, mutations, expressions):
		normalized_symbol = (symbol or "").strip()
		if gene and gene.symbol:
			normalized_symbol = gene.symbol

		hgnc_id = (gene.hgnc_id or "").replace("HGNC:", "") if gene else ""
		entrez_id = self._normalize_entrez_id(gene.entrez_id if gene else "")
		ensembl_id = self._normalize_ensembl_gene_id(gene.ensembl_gene_id if gene else "")

		if not entrez_id and mutations:
			for row in mutations:
				entrez_id = self._normalize_entrez_id(getattr(row, "entrez_gene_id", ""))
				if entrez_id:
					break

		if not ensembl_id and expressions:
			ensembl_id = self._extract_ensembl_from_expression(expressions[0])

		return {
			"hgnc": f"https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{hgnc_id}" if hgnc_id else "",
			"ensembl": f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={ensembl_id}" if ensembl_id else "",
			"ncbi": f"https://www.ncbi.nlm.nih.gov/gene/{entrez_id}" if entrez_id else "",
			"symbol_search": f"https://www.genecards.org/Search/Keyword?queryString={normalized_symbol}" if normalized_symbol else "",
		}

	def _enrich_rows_with_links(self, mutations, expressions):
		for row in mutations:
			file_id = getattr(getattr(row, "tcga_file", None), "file_id", "")
			entrez_id = self._normalize_entrez_id(getattr(row, "entrez_gene_id", ""))
			setattr(
				row,
				"gdc_file_url",
				f"https://portal.gdc.cancer.gov/files/{file_id}" if file_id else "",
			)
			setattr(
				row,
				"ncbi_gene_url",
				f"https://www.ncbi.nlm.nih.gov/gene/{entrez_id}" if entrez_id else "",
			)

		for row in expressions:
			file_id = getattr(getattr(row, "tcga_file", None), "file_id", "")
			ensembl_id = self._extract_ensembl_from_expression(row)
			setattr(
				row,
				"gdc_file_url",
				f"https://portal.gdc.cancer.gov/files/{file_id}" if file_id else "",
			)
			setattr(
				row,
				"ensembl_gene_url",
				f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={ensembl_id}" if ensembl_id else "",
			)

	def _parse_float(self, value):
		value = (value or "").strip()
		if not value:
			return None
		try:
			return float(value)
		except (TypeError, ValueError):
			return None

	def _apply_limit(self, mutations_qs, expressions_qs, row_limit):
		if row_limit == "all":
			return mutations_qs, expressions_qs
		if row_limit not in {"100", "1000"}:
			row_limit = "100"
		limit_value = int(row_limit)
		return mutations_qs[:limit_value], expressions_qs[:limit_value]

	def _build_raw_entries(self, mutations_qs, expressions_qs):
		file_values = {}
		for file_id, file_name, local_path in mutations_qs.values_list(
			"tcga_file__file_id", "tcga_file__file_name", "tcga_file__local_path"
		):
			file_values[file_id] = (file_name, local_path)
		for file_id, file_name, local_path in expressions_qs.values_list(
			"tcga_file__file_id", "tcga_file__file_name", "tcga_file__local_path"
		):
			file_values[file_id] = (file_name, local_path)

		raw_entries = []
		missing_raw_paths = []
		for file_id, (file_name, local_path) in sorted(
			file_values.items(),
			key=lambda item: (
				str(item[0] or ""),
				str(item[1][0] or ""),
			),
		):
			path_value = (local_path or "").strip()
			if not path_value:
				missing_raw_paths.append(f"{file_name} (sin local_path)")
				continue
			if not os.path.exists(path_value):
				missing_raw_paths.append(path_value)
				continue
			raw_entries.append(
				{
					"arcname": f"raw/{file_id}_{file_name}",
					"path": path_value,
				}
			)
		return raw_entries, missing_raw_paths

	def _export_metadata(self, symbol, exported_at, missing_raw_paths):
		lines = [
			f"symbol={symbol}",
			f"exported_at_utc={exported_at}",
			f"missing_raw_files={len(missing_raw_paths)}",
		]
		if missing_raw_paths:
			lines.append("missing_raw_paths=")
			lines.extend(missing_raw_paths)
		return "\n".join(lines) + "\n"

	def _mutation_rows_for_export(self, mutations_qs, symbol, exported_at):
		rows = []
		for row in mutations_qs.values(*self.MUTATION_EXPORT_FIELDS):
			row["export_symbol"] = symbol
			row["exported_at_utc"] = exported_at
			rows.append(row)
		return rows

	def _expression_rows_for_export(self, expressions_qs, symbol, exported_at):
		rows = []
		for row in expressions_qs.values(*self.EXPRESSION_EXPORT_FIELDS):
			row["export_symbol"] = symbol
			row["exported_at_utc"] = exported_at
			rows.append(row)
		return rows

	def _csv_response(self, filename, rows, columns):
		output = io.StringIO()
		writer = csv.DictWriter(output, fieldnames=columns)
		writer.writeheader()
		for row in rows:
			writer.writerow(row)
		response = HttpResponse(output.getvalue(), content_type="text/csv")
		response["Content-Disposition"] = f'attachment; filename="{filename}"'
		return response

	def _xlsx_response(self, filename, rows):
		df = pd.DataFrame(rows)
		output = io.BytesIO()
		with pd.ExcelWriter(output, engine="openpyxl") as writer:
			df.to_excel(writer, index=False, sheet_name="data")
		output.seek(0)
		response = HttpResponse(
			output.read(),
			content_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
		)
		response["Content-Disposition"] = f'attachment; filename="{filename}"'
		return response

	def _build_export_response(
		self,
		export_kind,
		gene,
		symbol,
		exported_at,
		mutations_qs,
		expressions_qs,
		raw_entries,
		missing_raw_paths,
		request,
	):
		mutation_rows = self._mutation_rows_for_export(mutations_qs, symbol, exported_at)
		expression_rows = self._expression_rows_for_export(expressions_qs, symbol, exported_at)
		base_name = f"{symbol}_{exported_at}"
		metadata_text = self._export_metadata(symbol, exported_at, missing_raw_paths)
		username = request.user.get_username() if request.user.is_authenticated else "anonymous"

		if export_kind == "mut_csv":
			logger.info("Gene export user=%s type=mut_csv symbol=%s rows=%d", username, symbol, len(mutation_rows))
			return self._csv_response(
				f"{base_name}_mutations.csv",
				mutation_rows,
				self.MUTATION_EXPORT_FIELDS + ["export_symbol", "exported_at_utc"],
			)

		if export_kind == "expr_csv":
			logger.info("Gene export user=%s type=expr_csv symbol=%s rows=%d", username, symbol, len(expression_rows))
			return self._csv_response(
				f"{base_name}_expression.csv",
				expression_rows,
				self.EXPRESSION_EXPORT_FIELDS + ["export_symbol", "exported_at_utc"],
			)

		if export_kind == "mut_xlsx":
			logger.info("Gene export user=%s type=mut_xlsx symbol=%s rows=%d", username, symbol, len(mutation_rows))
			return self._xlsx_response(f"{base_name}_mutations.xlsx", mutation_rows)

		if export_kind == "expr_xlsx":
			logger.info("Gene export user=%s type=expr_xlsx symbol=%s rows=%d", username, symbol, len(expression_rows))
			return self._xlsx_response(f"{base_name}_expression.xlsx", expression_rows)

		if export_kind == "raw_zip":
			logger.info("Gene export user=%s type=raw_zip symbol=%s files=%d", username, symbol, len(raw_entries))
			return self._raw_zip_response(f"{base_name}_raw.zip", raw_entries, metadata_text)

		if export_kind == "all_zip":
			logger.info(
				"Gene export user=%s type=all_zip symbol=%s mut_rows=%d expr_rows=%d raw_files=%d",
				username,
				symbol,
				len(mutation_rows),
				len(expression_rows),
				len(raw_entries),
			)
			return self._full_zip_response(
				filename=f"{base_name}_bundle.zip",
				mutation_rows=mutation_rows,
				expression_rows=expression_rows,
				raw_entries=raw_entries,
				metadata_text=metadata_text,
			)

		return None

	def _raw_zip_response(self, filename, raw_entries, metadata_text):
		output = io.BytesIO()
		with zipfile.ZipFile(output, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
			zf.writestr("metadata.txt", metadata_text)
			for entry in raw_entries:
				zf.write(entry["path"], arcname=entry["arcname"])
		output.seek(0)
		response = HttpResponse(output.read(), content_type="application/zip")
		response["Content-Disposition"] = f'attachment; filename="{filename}"'
		return response

	def _full_zip_response(self, filename, mutation_rows, expression_rows, raw_entries, metadata_text):
		output = io.BytesIO()
		with zipfile.ZipFile(output, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
			zf.writestr("metadata.txt", metadata_text)

			mutation_buffer = io.StringIO()
			mutation_writer = csv.DictWriter(
				mutation_buffer,
				fieldnames=self.MUTATION_EXPORT_FIELDS + ["export_symbol", "exported_at_utc"],
			)
			mutation_writer.writeheader()
			for row in mutation_rows:
				mutation_writer.writerow(row)
			zf.writestr("processed/mutations.csv", mutation_buffer.getvalue())

			expression_buffer = io.StringIO()
			expression_writer = csv.DictWriter(
				expression_buffer,
				fieldnames=self.EXPRESSION_EXPORT_FIELDS + ["export_symbol", "exported_at_utc"],
			)
			expression_writer.writeheader()
			for row in expression_rows:
				expression_writer.writerow(row)
			zf.writestr("processed/expression.csv", expression_buffer.getvalue())

			for entry in raw_entries:
				zf.write(entry["path"], arcname=entry["arcname"])

		output.seek(0)
		response = HttpResponse(output.read(), content_type="application/zip")
		response["Content-Disposition"] = f'attachment; filename="{filename}"'
		return response
