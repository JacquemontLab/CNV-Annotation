#!/usr/bin/env python
import polars as pl
import sys

input = sys.argv[1]
output = sys.argv[2]

# Expecting existing columns
# SampleID      Chr     Start   End     Type

# Convert column names to lowercase for case-insensitive mapping
with open(input) as f:
    header = f.readline().strip().split("\t")
col_map = {name.lower(): name for name in header}

# Scan csv and cast to schema
df  = pl.scan_csv(input,
                  separator="\t",
                  infer_schema_length=1000000)

# Create CNV_ID column
df = df.with_columns(
    pl.concat_str([
        pl.col(col_map["chr"]), 
        pl.col(col_map["start"]), 
        pl.col(col_map["end"]), 
        pl.col(col_map["type"])], separator="_").alias("CNV_ID")
)

# Defining output order, putting IDs in front
order = (["CNV_ID", "SampleID"] +
         [col for col in df.columns if col not in ["CNV_ID", "SampleID"]])
df = df.select(order)

df.sink_parquet(output)
