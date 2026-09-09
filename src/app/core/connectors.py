import os
import json
from typing import TypedDict

import requests
import pandas as pd
import GEOparse
from rich.console import Console
from rich.table import Table

console = Console(no_color=True, width=80)


class GEOStudyInfo(TypedDict):
    accession: str
    title: str
    n_samples: int
    taxon: str
    gdstype: str
    pdat: str


class GEOMetadataInfo(TypedDict):
    accession: str
    title: str
    platform: str
    n_samples: int
    samples_df: pd.DataFrame


class GEOConnector:
    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def search(self, text: str, num_results: int = 20) -> list[GEOStudyInfo]:
        # esearch devuelve ids internos de ncbi (no son los gse)
        params = {"db": "gds", "term": text, "retmax": num_results, "retmode": "json"}
        r = requests.get(f"{self.BASE_URL}/esearch.fcgi", params=params, timeout=15)
        r.raise_for_status()
        ids = r.json()["esearchresult"]["idlist"]

        if not ids:
            return []

        # summary de cada id
        params = {"db": "gds", "id": ",".join(ids), "retmode": "json"}
        r = requests.get(f"{self.BASE_URL}/esummary.fcgi", params=params, timeout=15)
        r.raise_for_status()
        result = r.json()["result"]

        studies = []
        for uid in ids:
            entry = result.get(uid, {})
            if entry.get("entrytype") != "GSE":  # solo series completas
                continue
            studies.append(GEOStudyInfo(
                accession=entry.get("accession", ""),
                title=entry.get("title", ""),
                n_samples=entry.get("n_samples", ""),
                taxon=entry.get("taxon", ""),
                gdstype=entry.get("gdstype", ""),
                pdat=entry.get("pdat", ""),
            ))

        return studies

    def get_series_metadata(self, gse_id: str) -> GEOMetadataInfo:
        # descarga el soft file de la serie; gsm.metadata tiene todos los campos
        console.log(f"Descargando metadatos de {gse_id} ...")
        gse = GEOparse.get_GEO(geo=gse_id, destdir="/tmp", silent=True)

        samples_info = []
        for gsm_id, gsm in gse.gsms.items():
            row = {"gsm": gsm_id}
            row.update(gsm.metadata)
            samples_info.append(row)

        return GEOMetadataInfo(
            accession=gse_id,
            title=gse.metadata.get("title", [""])[0],
            platform=gse.metadata.get("platform_id", [""])[0],
            n_samples=len(gse.gsms),
            samples_df=pd.json_normalize(samples_info),
        )


class TCGAConnector:
    BASE_URL = "https://api.gdc.cancer.gov"

    def list_projects(self, keyword: str = "breast") -> pd.DataFrame:
        # la api de projects no filtra bien por disease_type libre, asi que
        # descarga todos y filtra en pandas
        url = f"{self.BASE_URL}/projects"
        params = {
            "fields": "project_id,name,disease_type,primary_site,summary.case_count,summary.file_count",
            "format": "json",
            "size": 100,
        }
        r = requests.get(url, params=params, timeout=15)
        r.raise_for_status()
        hits = r.json()["data"]["hits"]
        df = pd.json_normalize(hits)
        if keyword and not df.empty:
            # solo en project_id y name para evitar falsos positivos
            id_match = df["project_id"].str.contains(keyword, case=False, na=False)
            name_match = df["name"].str.contains(keyword, case=False, na=False)
            df = df[id_match | name_match]
        return df.reset_index(drop=True)

    def list_files(self, project_id: str = "TCGA-BRCA",
                   data_type: str = "Gene Expression Quantification",
                   num_results: int = 10) -> pd.DataFrame:
        # descarga el fichero con GET /data/{file_id}
        url = f"{self.BASE_URL}/files"
        filters = {
            "op": "and",
            "content": [
                {"op": "=", "content": {"field": "cases.project.project_id", "value": project_id}},
                {"op": "=", "content": {"field": "data_type", "value": data_type}},
            ],
        }
        params = {
            "filters": json.dumps(filters),
            "fields": "file_id,file_name,data_type,experimental_strategy,file_size,cases.case_id",
            "format": "json",
            "size": num_results,
        }
        r = requests.get(url, params=params, timeout=15)
        r.raise_for_status()
        hits = r.json()["data"]["hits"]
        return pd.json_normalize(hits)

    def download_file(
        self,
        file_id: str,
        dest_folder: str = "./tcga_downloads",
        timeout: int = 300,
    ) -> dict:
        # descarga un fichero tcga por file_id via GET /data/{file_id}
        url = f"{self.BASE_URL}/data/{file_id}"
        os.makedirs(dest_folder, exist_ok=True)

        r = requests.get(url, timeout=timeout, stream=True)
        r.raise_for_status()

        filename = self._extract_filename_from_response(r, fallback=file_id)
        filepath = os.path.join(dest_folder, filename)
        downloaded_bytes = 0

        with open(filepath, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                if not chunk:
                    continue
                f.write(chunk)
                downloaded_bytes += len(chunk)

        console.log(f"Archivo descargado: {filepath}")
        return {
            "file_id": file_id,
            "file_name": filename,
            "local_path": filepath,
            "downloaded_bytes": downloaded_bytes,
        }

    def download_files(
        self,
        file_list: list[str],
        dest_folder: str = "./tcga_downloads",
        related_files: bool = False,
        timeout: int = 300,
    ) -> list[dict]:
        # descarga una lista de file_ids, uno a uno
        if related_files:
            console.log("related_files no se usa en descarga por file_id; se ignora.")

        file_list = [file_list] if isinstance(file_list, str) else file_list
        results = []
        for file_id in file_list:
            results.append(self.download_file(file_id=file_id, dest_folder=dest_folder, timeout=timeout))
        return results

    def _extract_filename_from_response(self, response: requests.Response, fallback: str = "download") -> str:
        content_disposition = response.headers.get("Content-Disposition", "")
        if "filename=" in content_disposition:
            return content_disposition.split("filename=")[-1].strip('"')
        return fallback

    def get_clinical_summary(self, project_id: str = "TCGA-BRCA",
                              num_results: int = 20) -> pd.DataFrame:
        # datos clinicos mas completos estan en los ficheros de clinical supplement
        url = f"{self.BASE_URL}/cases"
        filters = {
            "op": "=",
            "content": {"field": "project.project_id", "value": project_id},
        }
        params = {
            "filters": json.dumps(filters),
            "fields": "case_id,primary_site,disease_type,demographic.gender,demographic.age_at_index",
            "format": "json",
            "size": num_results,
        }
        r = requests.get(url, params=params, timeout=15)
        r.raise_for_status()
        hits = r.json()["data"]["hits"]
        return pd.json_normalize(hits)


class HGNCConnector:
    BASE_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"

    def download_complete_set(self, dest_path: str = "/tmp/hgnc_complete_set.txt") -> pd.DataFrame:
        # descarga el tsv completo de hgnc y lo devuelve como dataframe
        if not os.path.exists(dest_path):
            console.log(f"Descargando HGNC complete set en {dest_path}...")
            r = requests.get(self.BASE_URL, timeout=60, stream=True)
            r.raise_for_status()
            with open(dest_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
            console.log("Descarga completada.")
        else:
            console.log(f"Usando fichero HGNC existente: {dest_path}")

        df = pd.read_csv(dest_path, sep="\t", dtype=str, low_memory=False)
        df = df.where(df != "", other=None)
        return df


def df_to_rich_table(df: pd.DataFrame, title: str) -> Table:
    # convierte un dataframe en una rich table en blanco y negro
    table = Table(title=title, header_style="bold", show_lines=True)
    for col in df.columns:
        table.add_column(col, overflow="fold")
    for _, row in df.iterrows():
        table.add_row(*[str(v) for v in row])
    return table


# main de pruebas
if __name__ == "__main__":
    geo = GEOConnector()
    studies = geo.search("breast cancer[title] AND expression profiling by array[DataSet Type]", num_results=10)
    df_geo = pd.DataFrame(studies)[["accession", "n_samples", "pdat", "title"]]
    console.print(df_to_rich_table(df_geo, "Series GEO - cancer de mama"))

    tcga = TCGAConnector()
    df_proj = tcga.list_projects(keyword="breast")
    show_cols = [c for c in ["project_id", "name", "summary.case_count", "summary.file_count"] if c in df_proj.columns]
    console.print(df_to_rich_table(df_proj[show_cols], "Proyectos TCGA - mama (breast)"))

    df_files = tcga.list_files(num_results=5)
    tcga.download_files(df_files["file_id"].tolist()[0], dest_folder="./downloads", related_files=True)
    show_fcols = [c for c in ["file_name", "data_type", "experimental_strategy", "file_size"] if c in df_files.columns]
    console.print(df_to_rich_table(df_files[show_fcols], "Ficheros TCGA-BRCA - Expresion genetica"))

    df_clin = tcga.get_clinical_summary(num_results=20)
    console.print(df_to_rich_table(df_clin, "Casos TCGA-BRCA - datos clinicos"))
