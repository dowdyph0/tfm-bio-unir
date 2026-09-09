# models.py
from django.db import models


class HGNCGene(models.Model):
    # tabla de referencia de genes humanos, cargada desde el tsv completo de hgnc
    hgnc_id = models.CharField(max_length=20, unique=True)
    symbol = models.CharField(max_length=50, db_index=True)
    name = models.TextField(blank=True)
    locus_group = models.CharField(max_length=100, blank=True)
    locus_type = models.CharField(max_length=100, blank=True)
    status = models.CharField(max_length=20, blank=True)
    location = models.CharField(max_length=50, blank=True)

    # ids cruzados para mapear geo y tcga
    entrez_id = models.CharField(max_length=20, blank=True, db_index=True)
    ensembl_gene_id = models.CharField(max_length=20, blank=True, db_index=True)
    refseq_accession = models.CharField(max_length=50, blank=True)
    uniprot_ids = models.CharField(max_length=100, blank=True)

    # simbolos alternativos
    alias_symbol = models.TextField(blank=True)
    prev_symbol = models.TextField(blank=True)

    class Meta:
        db_table = "hgnc_gene"
        indexes = [
            models.Index(fields=["symbol"]),
            models.Index(fields=["ensembl_gene_id"]),
            models.Index(fields=["entrez_id"]),
        ]

    def __str__(self):
        return f"{self.symbol} ({self.hgnc_id})"


class GEOSeries(models.Model):
    accession = models.CharField(max_length=20, unique=True)
    title = models.TextField()
    summary = models.TextField(blank=True)
    gds_type = models.CharField(max_length=100, blank=True)
    taxon = models.CharField(max_length=100, blank=True)
    platform_id = models.CharField(max_length=20, blank=True)
    n_samples = models.IntegerField(default=0)
    publication_date = models.CharField(max_length=10, blank=True)
    supplementary = models.CharField(max_length=200, blank=True)
    ftp_link = models.URLField(max_length=300, blank=True)
    bioproject = models.CharField(max_length=30, blank=True)

    platform_title = models.CharField(max_length=200, blank=True)

    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        db_table = "geo_series"

    def __str__(self):
        return self.accession


class GEOSample(models.Model):
    series = models.ForeignKey(GEOSeries, on_delete=models.CASCADE,
                               related_name="samples")
    gsm_id = models.CharField(max_length=20, unique=True)
    title = models.CharField(max_length=300, blank=True)
    source_name = models.CharField(max_length=200, blank=True)
    organism = models.CharField(max_length=100, blank=True)
    platform_id = models.CharField(max_length=20, blank=True)

    # caracteristicas clinicas heterogeneas en json
    clinical_metadata = models.JSONField(default=dict, blank=True)

    class Meta:
        db_table = "geo_sample"

    def __str__(self):
        return self.gsm_id


class GEOExpressionRecord(models.Model):
    # valor de expresion de una muestra geo para un gen hgnc
    geo_sample = models.ForeignKey(GEOSample, on_delete=models.CASCADE,
                                   related_name="expression_records")
    hgnc_gene = models.ForeignKey(HGNCGene, on_delete=models.CASCADE,
                                  related_name="geo_expression_records",
                                  null=True, blank=True)
    gene_symbol = models.CharField(max_length=50, db_index=True)
    value = models.FloatField()
    zscore = models.FloatField(null=True, blank=True)

    class Meta:
        db_table = "geo_expression_record"
        indexes = [
            models.Index(fields=["gene_symbol"]),
            models.Index(fields=["geo_sample", "gene_symbol"]),
        ]
        constraints = [
            models.UniqueConstraint(
                fields=["geo_sample", "gene_symbol"],
                name="uq_geo_expr_sample_gene"
            )
        ]

    def __str__(self):
        return f"{self.geo_sample_id}:{self.gene_symbol}={self.value:.3f}"


class TCGAProject(models.Model):
    project_id = models.CharField(max_length=50, unique=True)
    name = models.TextField()
    # disease_type y primary_site son listas en la api
    disease_type = models.JSONField(default=list)
    primary_site = models.JSONField(default=list)
    case_count = models.IntegerField(default=0)
    file_count = models.IntegerField(default=0)

    class Meta:
        db_table = "tcga_project"

    def __str__(self):
        return self.project_id


class TCGACase(models.Model):
    project = models.ForeignKey(TCGAProject, on_delete=models.CASCADE,
                                related_name="cases",
                                null=True, blank=True,
                                to_field="project_id")
    case_id = models.CharField(max_length=50, unique=True)
    primary_site = models.CharField(max_length=100, blank=True)
    disease_type = models.CharField(max_length=200, blank=True)
    gender = models.CharField(max_length=20, blank=True)
    age_at_index = models.IntegerField(null=True, blank=True)
    molecular_subtype = models.CharField(max_length=50, blank=True)  # pam50 real
    predicted_subtype = models.CharField(max_length=50, blank=True)  # pam50 predicho

    class Meta:
        db_table = "tcga_case"

    def __str__(self):
        return self.case_id


class TCGAFile(models.Model):
    cases = models.ManyToManyField(TCGACase, related_name="files", blank=True)
    file_id = models.CharField(max_length=50, unique=True)
    file_name = models.CharField(max_length=300)
    data_type = models.CharField(max_length=100, blank=True)
    experimental_strategy = models.CharField(max_length=100, blank=True)
    file_size = models.BigIntegerField(default=0)
    local_path = models.CharField(max_length=500, blank=True)

    class Meta:
        db_table = "tcga_file"

    def __str__(self):
        return self.file_name


class DownloadStatus(models.Model):
    # estado de descarga/procesado de un tcga file
    STATUS_PENDING = "pending"
    STATUS_DOWNLOADING = "downloading"
    STATUS_DOWNLOADED = "downloaded"
    STATUS_PARSED = "parsed"
    STATUS_FAILED = "failed"

    STATUS_CHOICES = [
        (STATUS_PENDING, "Pending"),
        (STATUS_DOWNLOADING, "Downloading"),
        (STATUS_DOWNLOADED, "Downloaded"),
        (STATUS_PARSED, "Parsed"),
        (STATUS_FAILED, "Failed"),
    ]

    tcga_file = models.OneToOneField(
        TCGAFile,
        on_delete=models.CASCADE,
        related_name="download_status",
    )
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default=STATUS_PENDING)
    started_at = models.DateTimeField(null=True, blank=True)
    finished_at = models.DateTimeField(null=True, blank=True)
    downloaded_bytes = models.BigIntegerField(default=0)
    local_path = models.CharField(max_length=500, blank=True)
    parser_used = models.CharField(max_length=100, blank=True)
    parse_summary = models.JSONField(default=dict, blank=True)
    error_message = models.TextField(blank=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        db_table = "download_status"
        indexes = [
            models.Index(fields=["status"]),
            models.Index(fields=["updated_at"]),
        ]

    def __str__(self):
        return f"{self.tcga_file.file_id}: {self.status}"


class TCGAMutation(models.Model):
    # registros parseados desde ficheros maf de tcga
    tcga_file = models.ForeignKey(TCGAFile, on_delete=models.CASCADE, related_name="mutations")
    hgnc_gene = models.ForeignKey(HGNCGene, null=True, blank=True, on_delete=models.SET_NULL, related_name="tcga_mutations")

    hugo_symbol = models.CharField(max_length=50, blank=True, db_index=True)
    entrez_gene_id = models.CharField(max_length=20, blank=True, db_index=True)
    variant_classification = models.CharField(max_length=80, blank=True, db_index=True)
    variant_type = models.CharField(max_length=40, blank=True)
    chromosome = models.CharField(max_length=20, blank=True)
    start_position = models.BigIntegerField(null=True, blank=True)
    end_position = models.BigIntegerField(null=True, blank=True)
    reference_allele = models.CharField(max_length=200, blank=True)
    tumor_seq_allele2 = models.CharField(max_length=200, blank=True)
    tumor_sample_barcode = models.CharField(max_length=80, blank=True, db_index=True)
    raw_row = models.JSONField(default=dict, blank=True)

    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "tcga_mutation"
        indexes = [
            models.Index(fields=["hugo_symbol"]),
            models.Index(fields=["entrez_gene_id"]),
            models.Index(fields=["variant_classification"]),
        ]


class TCGAExpressionRecord(models.Model):
    # registros parseados desde ficheros de expresion de tcga
    tcga_file = models.ForeignKey(TCGAFile, on_delete=models.CASCADE, related_name="expression_records")
    hgnc_gene = models.ForeignKey(HGNCGene, null=True, blank=True, on_delete=models.SET_NULL, related_name="tcga_expression_records")

    gene_id = models.CharField(max_length=40, blank=True, db_index=True)
    gene_id_versionless = models.CharField(max_length=30, blank=True, db_index=True)
    gene_name = models.CharField(max_length=120, blank=True, db_index=True)
    gene_type = models.CharField(max_length=120, blank=True)

    unstranded = models.FloatField(null=True, blank=True)
    stranded_first = models.FloatField(null=True, blank=True)
    stranded_second = models.FloatField(null=True, blank=True)
    tpm_unstranded = models.FloatField(null=True, blank=True)
    fpkm_unstranded = models.FloatField(null=True, blank=True)
    fpkm_uq_unstranded = models.FloatField(null=True, blank=True)
    raw_row = models.JSONField(default=dict, blank=True)

    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = "tcga_expression_record"
        indexes = [
            models.Index(fields=["gene_id"]),
            models.Index(fields=["gene_id_versionless"]),
            models.Index(fields=["gene_name"]),
        ]


class MLModelRun(models.Model):
    # registro de cada ejecucion del clasificador pam50
    run_date = models.DateTimeField(auto_now_add=True)
    algorithm = models.CharField(max_length=80)
    n_samples = models.IntegerField(default=0)
    n_features = models.IntegerField(default=0)
    test_size = models.FloatField(default=0.2)
    random_state = models.IntegerField(default=42)
    hyperparams = models.JSONField(default=dict, blank=True)
    accuracy = models.FloatField(null=True, blank=True)
    f1_macro = models.FloatField(null=True, blank=True)
    f1_weighted = models.FloatField(null=True, blank=True)
    confusion_matrix = models.JSONField(default=list, blank=True)
    class_labels = models.JSONField(default=list, blank=True)
    classification_report = models.TextField(blank=True)
    notes = models.TextField(blank=True)

    class Meta:
        db_table = "ml_model_run"
        ordering = ["-run_date"]

    def __str__(self):
        return f"{self.algorithm} {self.run_date:%Y-%m-%d %H:%M} acc={self.accuracy:.3f}"


class MLGeneImportance(models.Model):
    # importancia de cada gen en la clasificacion pam50
    run = models.ForeignKey(MLModelRun, on_delete=models.CASCADE,
                            related_name="gene_importances")
    hgnc_gene = models.ForeignKey(HGNCGene, null=True, blank=True,
                                  on_delete=models.SET_NULL)
    gene_name = models.CharField(max_length=120, db_index=True)
    gene_id = models.CharField(max_length=40, blank=True)
    importance = models.FloatField()
    rank = models.IntegerField()

    class Meta:
        db_table = "ml_gene_importance"
        ordering = ["rank"]
        indexes = [
            models.Index(fields=["run", "rank"]),
            models.Index(fields=["gene_name"]),
        ]

    def __str__(self):
        return f"#{self.rank} {self.gene_name} ({self.importance:.4f})"