#!/usr/bin/env python
"""
Script to preprocess a tab-delimited CNV file by extracting and renaming positional columns
(Chromosome, Start, End, CNV type), adding a placeholder Strand column, and exporting a
unique, sorted list of CNV regions in TSV format without a header.

Usage:
    python script.py <input_file.tsv> <output_file.tsv>

Arguments:
    input_file.tsv   Input CNV file with at least 5 columns (e.g., SampleID, Chr, Start, End, TYPE)
    output_file.tsv  Output TSV file with columns: Chr, Start, End, TYPE, Strand

Dependencies:
    - polars >= 0.20
"""


import polars as pl
import polars.selectors as cs
import sys


input_file = sys.argv[1]
output_file_bed = sys.argv[2]
output_file_parquet = sys.argv[3]


df = pl.scan_csv(input_file, separator = "\t", infer_schema_length=10000)

# Expecting existing columns
# SampleID  Chr     Start   End     Type    

# Convert column names to lowercase for case-insensitive mapping
with open(input_file) as f:
    header = f.readline().strip().split("\t")
col_map = {name.lower(): name for name in header}

# Get columns based on name, case-insensitive
df = df.select([
    pl.col(col_map["chr"]).alias("Chr"),
    pl.col(col_map["start"]).alias("Start"),
    pl.col(col_map["end"]).alias("End"),
    pl.col(col_map["type"]).alias("Type"),
    pl.lit(".").alias("Strand")
])


# --- Add CNV_ID (Chr_Start_End_Type) ---
df = df.with_columns(
    pl.concat_str([
        pl.col("Chr"),
        pl.col("Start").cast(pl.Utf8),
        pl.col("End").cast(pl.Utf8),
        pl.col("Type")
    ], separator="_").alias("CNV_ID")
)

# --- Keep only unique CNV_IDs and sort by Chr, Start, End ---
df_unique = df.select(["Chr", "Start", "End", "Type", "Strand", "CNV_ID"]) \
              .unique(subset="CNV_ID", keep="any") \
              .sort(["Chr", "Start", "End"])

# --- Write TSV without header ---
df_unique.select(["Chr", "Start", "End", "Type", "Strand"]) \
         .sink_csv(output_file_bed, separator="\t", include_header=False)

# --- Write Parquet (full data with CNV_ID) ---
df_unique.sink_parquet(output_file_parquet)