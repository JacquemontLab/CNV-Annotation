#!/usr/bin/env python
"""
merge_cnv_with_region.py

This script merges a CNV dataset with a genomic region overlap file.

Usage:
    python merge_cnv_with_region.py <cnv_file> <region_file> <output_parquet>

Arguments:
    cnv_file     : Path to the input CNV TSV file
    region_file  : Path to the region overlap TSV file
    output       : Path to the output Parquet file

"""

import polars as pl
import sys

cnv_file = sys.argv[1]
region_file = sys.argv[2]
output = sys.argv[3]

# --- Load CNV file ---
with open(cnv_file) as f:
    header = f.readline().strip().split("\t")
col_map = {name.lower(): name for name in header}

df = pl.read_csv(
    cnv_file,
    separator="\t",
    infer_schema_length=1000000,
    schema_overrides={"SampleID": pl.Utf8}  # Ensure SampleID stays string
)

# Build CNV_ID in same format as region_overlap
df = df.with_columns(
    pl.concat_str([
        pl.col(col_map["chr"]), 
        pl.col(col_map["start"]).cast(pl.Utf8), 
        pl.col(col_map["end"]).cast(pl.Utf8), 
        pl.col(col_map["type"])
    ], separator="_").alias("CNV_ID")
)

# --- Load overlap file ---
region_df = pl.read_csv(
    region_file,
    separator="\t",
    infer_schema_length=1000000
)

# --- Merge on CNV_ID ---
df = df.join(region_df, on="CNV_ID", how="left")

# --- Column order (IDs in front) ---
order = (["CNV_ID", "SampleID"] +
         [col for col in df.columns if col not in ["CNV_ID", "SampleID"]])
df = df.select(order)

# --- Save ---
df.write_parquet(output, compression="zstd")
