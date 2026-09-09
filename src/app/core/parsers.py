import os
from typing import Any

import pandas as pd


def parse_maf(file_path: str) -> tuple[pd.DataFrame, dict[str, Any]]:
    # lee un fichero maf (.maf o .maf.gz) y devuelve dataframe + resumen
    df = pd.read_csv(file_path, sep="\t", comment="#", low_memory=False)

    summary = {
        "rows": int(len(df)),
        "columns": list(df.columns),
        "unique_genes": int(df["Hugo_Symbol"].nunique()) if "Hugo_Symbol" in df.columns else 0,
        "unique_samples": int(df["Tumor_Sample_Barcode"].nunique()) if "Tumor_Sample_Barcode" in df.columns else 0,
    }
    return df, summary


def parse_expression_quantification(file_path: str) -> tuple[pd.DataFrame, dict[str, Any]]:
    # los tsv de tcga star counts incluyen metadatos en lineas que empiezan por '#'
    df = pd.read_csv(file_path, sep="\t", low_memory=False, comment="#")

    quant_columns = [
        c for c in [
            "unstranded",
            "stranded_first",
            "stranded_second",
            "tpm_unstranded",
            "fpkm_unstranded",
            "fpkm_uq_unstranded",
        ] if c in df.columns
    ]

    summary = {
        "rows": int(len(df)),
        "columns": list(df.columns),
        "n_quant_columns": len(quant_columns),
        "quant_columns": quant_columns,
        "n_gene_ids": int(df["gene_id"].nunique()) if "gene_id" in df.columns else 0,
    }
    return df, summary


def infer_parser(file_name: str, data_type: str = "") -> str:
    lowered = file_name.lower()
    data_type_lower = (data_type or "").lower()

    if lowered.endswith(".maf") or lowered.endswith(".maf.gz"):
        return "maf"
    if "masked_somatic_mutation" in lowered:
        return "maf"

    if "gene_expression_quantification" in data_type_lower:
        return "expression"
    if "expression quantification" in data_type_lower:
        return "expression"
    if lowered.endswith(".rna_seq.augmented_star_gene_counts.tsv"):
        return "expression"
    if lowered.endswith(".tsv") and "gene_counts" in lowered:
        return "expression"

    return "unknown"


def parse_tcga_file(
    file_path: str,
    file_name: str,
    data_type: str = "",
    include_dataframe: bool = False,
) -> tuple[str, dict[str, Any], pd.DataFrame | None]:
    # parsea un fichero tcga conocido y devuelve parser usado + resumen (+ dataframe opcional)
    parser_name = infer_parser(file_name=file_name, data_type=data_type)

    if parser_name == "maf":
        df, summary = parse_maf(file_path)
        return parser_name, summary, df if include_dataframe else None
    if parser_name == "expression":
        df, summary = parse_expression_quantification(file_path)
        return parser_name, summary, df if include_dataframe else None

    file_size = os.path.getsize(file_path) if os.path.exists(file_path) else 0
    return "none", {
        "rows": None,
        "columns": [],
        "note": "No parser available for this file type",
        "file_size_bytes": int(file_size),
    }, None
