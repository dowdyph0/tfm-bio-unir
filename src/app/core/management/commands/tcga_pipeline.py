from django.core.management.base import BaseCommand

from core.models import TCGAFile
from core.pipeline import TCGADownloadPipeline


class Command(BaseCommand):
    help = "Descarga y parsea ficheros TCGA con tracking en DownloadStatus"

    def add_arguments(self, parser):
        parser.add_argument("--project-id", default="TCGA-BRCA")
        parser.add_argument(
            "--data-types",
            nargs="+",
            default=[
                "Masked Somatic Mutation",
                "Gene Expression Quantification",
            ],
        )
        parser.add_argument("--num-results", type=int, default=2)
        parser.add_argument("--dest-folder", default="/app/downloads")
        parser.add_argument("--skip-parse", action="store_true")

    def handle(self, *args, **options):
        project_id = options["project_id"]
        data_types = options["data_types"]
        num_results = options["num_results"]
        dest_folder = options["dest_folder"]
        parse_downloaded = not options["skip_parse"]

        pipeline = TCGADownloadPipeline()

        enqueued_total = 0
        for data_type in data_types:
            n = pipeline.enqueue_project_files(
                project_id=project_id,
                data_type=data_type,
                num_results=num_results,
            )
            enqueued_total += n
            self.stdout.write(self.style.SUCCESS(f"[enqueue] {data_type}: {n} ficheros"))

        qs = TCGAFile.objects.filter(data_type__in=data_types).order_by("file_name")
        result = pipeline.run_for_existing_files(
            queryset=qs,
            dest_folder=dest_folder,
            parse_downloaded=parse_downloaded,
        )

        self.stdout.write(self.style.SUCCESS("Pipeline TCGA finalizado"))
        self.stdout.write(f"Enqueued: {enqueued_total}")
        self.stdout.write(f"Total: {result['total']}")
        self.stdout.write(f"Downloaded: {result['downloaded']}")
        self.stdout.write(f"Parsed: {result['parsed']}")
        self.stdout.write(f"Failed: {result['failed']}")
        self.stdout.write(f"Skipped existing: {result['skipped']}")
