#!/usr/bin/env python
import polars as pl
import sys

input = sys.argv[1]
output = sys.argv[2]

#Expects CNV file to at least contain the following columns in the proceeding order:
# SampleID      Chr     Start   End     Type

#scan csv and cast to schema
df  = pl.scan_csv(input,
                  separator="\t",
                  infer_schema_length=50000)


#Creating CNV_ID Column
df = df.with_columns(
        pl.concat_str([pl.col(col) for col in df.collect_schema().names()[1:5]],separator = "_").alias("CNV_ID")
)

#Defining output order, putting IDs in front
order = (["CNV_ID", "SampleID"] +
         [col for col in df.collect_schema().names() if col not in ["CNV_ID", "SampleID"]])

df = df.select(order)
df.sink_parquet(output)
