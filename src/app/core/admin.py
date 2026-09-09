from django.contrib import admin

from django.contrib.admin import SimpleListFilter

from .pipeline import TCGADownloadPipeline
from .models import (
    HGNCGene,
    GEOSeries,
    GEOSample,
    TCGAProject,
    TCGACase,
    TCGAFile,
    DownloadStatus,
    TCGAMutation,
    TCGAExpressionRecord,
    MLModelRun,
    MLGeneImportance,
)


class ProjectIdFilter(SimpleListFilter):
    title = "Project"
    parameter_name = "project_id"

    def lookups(self, request, model_admin):
        projects = (
            TCGAProject.objects.order_by("project_id")
            .values_list("project_id", flat=True)
            .distinct()
        )
        return [(p, p) for p in projects if p]

    def queryset(self, request, queryset):
        value = self.value()
        if not value:
            return queryset
        return queryset.filter(tcga_file__cases__project__project_id=value).distinct()





@admin.register(HGNCGene)
class HGNCGeneAdmin(admin.ModelAdmin):
    list_display = ("symbol", "hgnc_id", "name", "locus_type", "entrez_id", "ensembl_gene_id")
    list_filter = ("locus_group", "locus_type", "status")
    search_fields = ("symbol", "hgnc_id", "name", "alias_symbol", "entrez_id", "ensembl_gene_id")
    ordering = ("symbol",)
    readonly_fields = ("hgnc_id",)


@admin.register(GEOSeries)
class GEOSeriesAdmin(admin.ModelAdmin):
    list_display = ("accession", "title", "taxon", "gds_type", "n_samples", "publication_date")
    list_filter = ("taxon", "gds_type")
    search_fields = ("accession", "title", "summary", "bio_project")
    ordering = ("-publication_date",)
    readonly_fields = ("accession", "created_at", "updated_at")


@admin.register(GEOSample)
class GEOSampleAdmin(admin.ModelAdmin):
    list_display = ("gsm_id", "title", "organism", "platform_id", "series")
    list_filter = ("organism", "platform_id")
    search_fields = ("gsm_id", "title", "source_name")
    ordering = ("gsm_id",)
    readonly_fields = ("gsm_id",)
    raw_id_fields = ("series",)  # evita cargar todo el dropdown si hay muchas series


@admin.register(TCGAProject)
class TCGAProjectAdmin(admin.ModelAdmin):
    list_display = ("project_id", "name", "case_count", "file_count")
    search_fields = ("project_id", "name")
    ordering = ("project_id",)
    readonly_fields = ("project_id",)


@admin.register(TCGACase)
class TCGACaseAdmin(admin.ModelAdmin):
    list_display = ("case_id", "project", "primary_site", "disease_type", "gender", "age_at_index")
    list_filter = ("primary_site", "gender")
    search_fields = ("case_id", "primary_site", "disease_type")
    ordering = ("case_id",)
    readonly_fields = ("case_id",)
    raw_id_fields = ("project",)


@admin.register(TCGAFile)
class TCGAFileAdmin(admin.ModelAdmin):
    list_display = ("file_id", "file_name", "data_type", "experimental_strategy", "file_size_mb")
    list_filter = ("data_type", "experimental_strategy")
    search_fields = ("file_id", "file_name")
    ordering = ("file_name",)
    readonly_fields = ("file_id",)
    filter_horizontal = ("cases",)  # widget mas usable para el m2m

    @admin.display(description="Size (MB)")
    def file_size_mb(self, obj):
        if obj.file_size:
            return f"{obj.file_size / 1_048_576:.2f} MB"
        return "-"


@admin.register(DownloadStatus)
class DownloadStatusAdmin(admin.ModelAdmin):
    list_display = (
        "tcga_file",
        "status",
        "downloaded_bytes_mb",
        "parser_used",
        "updated_at",
    )
    list_filter = ("status", "parser_used")
    search_fields = (
        "tcga_file__file_id",
        "tcga_file__file_name",
        "local_path",
    )
    ordering = ("-updated_at",)
    raw_id_fields = ("tcga_file",)
    list_select_related = ("tcga_file",)
    list_per_page = 50
    show_full_result_count = False
    actions = ("reparse_selected", "redownload_and_parse_selected", "mark_pending")
    readonly_fields = (
        "tcga_file",
        "status",
        "started_at",
        "finished_at",
        "downloaded_bytes",
        "local_path",
        "parser_used",
        "parse_summary",
        "error_message",
        "created_at",
        "updated_at",
    )

    @admin.action(description="Re-parse selected files")
    def reparse_selected(self, request, queryset):
        pipeline = TCGADownloadPipeline()
        ok = 0
        failed = 0
        for status_obj in queryset.select_related("tcga_file"):
            if pipeline._parse_file(status_obj.tcga_file, status_obj):
                ok += 1
            else:
                failed += 1
        self.message_user(request, f"Re-parse finished: ok={ok}, failed={failed}")

    @admin.action(description="Re-download + parse selected files")
    def redownload_and_parse_selected(self, request, queryset):
        pipeline = TCGADownloadPipeline()
        file_qs = TCGAFile.objects.filter(id__in=queryset.values_list("tcga_file_id", flat=True))
        summary = pipeline.run_for_existing_files(
            queryset=file_qs,
            dest_folder="/app/downloads",
            parse_downloaded=True,
        )
        self.message_user(
            request,
            "Pipeline run: "
            f"total={summary['total']} parsed={summary['parsed']} failed={summary['failed']}",
        )

    @admin.action(description="Mark selected as pending")
    def mark_pending(self, request, queryset):
        updated = queryset.update(
            status=DownloadStatus.STATUS_PENDING,
            error_message="",
            parser_used="",
            parse_summary={},
        )
        self.message_user(request, f"Updated {updated} records to pending")

    @admin.display(description="Downloaded (MB)")
    def downloaded_bytes_mb(self, obj):
        if obj.downloaded_bytes:
            return f"{obj.downloaded_bytes / 1_048_576:.2f} MB"
        return "-"


@admin.register(TCGAMutation)
class TCGAMutationAdmin(admin.ModelAdmin):
    list_display = (
        "tcga_file",
        "hugo_symbol",
        "variant_classification",
        "chromosome",
        "start_position",
    )
    list_filter = (ProjectIdFilter, "variant_classification", "chromosome")
    search_fields = (
        "tcga_file__file_name",
        "hugo_symbol",
        "entrez_gene_id",
        "tumor_sample_barcode",
    )
    ordering = ("tcga_file", "hugo_symbol")
    raw_id_fields = ("tcga_file", "hgnc_gene")
    list_select_related = ("tcga_file", "hgnc_gene")
    list_per_page = 100
    show_full_result_count = False
    readonly_fields = (
        "tcga_file",
        "hgnc_gene",
        "hugo_symbol",
        "entrez_gene_id",
        "variant_classification",
        "variant_type",
        "chromosome",
        "start_position",
        "end_position",
        "reference_allele",
        "tumor_seq_allele2",
        "tumor_sample_barcode",
        "raw_row",
        "created_at",
    )

    def has_add_permission(self, request):
        return False


@admin.register(TCGAExpressionRecord)
class TCGAExpressionRecordAdmin(admin.ModelAdmin):
    list_display = (
        "tcga_file",
        "gene_id",
        "gene_name",
        "tpm_unstranded",
        "fpkm_uq_unstranded",
    )
    list_filter = (ProjectIdFilter, "gene_type")
    search_fields = (
        "tcga_file__file_name",
        "gene_id",
        "gene_name",
    )
    ordering = ("tcga_file", "gene_id")
    raw_id_fields = ("tcga_file", "hgnc_gene")
    list_select_related = ("tcga_file", "hgnc_gene")
    list_per_page = 100
    show_full_result_count = False
    readonly_fields = (
        "tcga_file",
        "hgnc_gene",
        "gene_id",
        "gene_id_versionless",
        "gene_name",
        "gene_type",
        "unstranded",
        "stranded_first",
        "stranded_second",
        "tpm_unstranded",
        "fpkm_unstranded",
        "fpkm_uq_unstranded",
        "raw_row",
        "created_at",
    )

    def has_add_permission(self, request):
        return False

class MLGeneImportanceInline(admin.TabularInline):
    model = MLGeneImportance
    fields = ("rank", "gene_name", "gene_id", "importance", "hgnc_gene")
    readonly_fields = ("rank", "gene_name", "gene_id", "importance", "hgnc_gene")
    extra = 0
    max_num = 50
    ordering = ("rank",)
    can_delete = False

    def has_add_permission(self, request, obj=None):
        return False


@admin.register(MLModelRun)
class MLModelRunAdmin(admin.ModelAdmin):
    list_display = ("run_date", "algorithm", "n_samples", "n_features", "accuracy", "f1_macro", "f1_weighted")
    readonly_fields = ("run_date", "algorithm", "n_samples", "n_features", "test_size", "random_state", "hyperparams", "accuracy", "f1_macro", "f1_weighted", "confusion_matrix", "class_labels", "classification_report")
    inlines = [MLGeneImportanceInline]
    ordering = ("-run_date",)

    def has_add_permission(self, request):
        return False


@admin.register(MLGeneImportance)
class MLGeneImportanceAdmin(admin.ModelAdmin):
    list_display = ("rank", "gene_name", "gene_id", "importance", "run")
    list_filter = ("run",)
    search_fields = ("gene_name", "gene_id")
    ordering = ("run", "rank")
    raw_id_fields = ("run", "hgnc_gene")
    list_per_page = 100
    show_full_result_count = False
    readonly_fields = ("run", "gene_name", "gene_id", "importance", "rank", "hgnc_gene")

    def has_add_permission(self, request):
        return False
