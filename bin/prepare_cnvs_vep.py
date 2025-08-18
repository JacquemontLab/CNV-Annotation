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
output_file = sys.argv[2]


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

# Final selection, deduplication, sorting and output
(df.select(["Chr", "Start", "End", "Type", "Strand"])
   .unique(keep="any")
   .sort(by=["Chr", "Start", "End"])
   .sink_csv(output_file, separator="\t", include_header=False))