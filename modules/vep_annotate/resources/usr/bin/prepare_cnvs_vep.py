#!/usr/bin/env python

import polars as pl
import polars.selectors as cs
import sys


input_file = sys.argv[1]
output_file = sys.argv[2]


df = pl.scan_csv(input_file, separator = "\t", infer_schema_length=10000)


#positional selection, expects
# SampleID  Chr     Start   End     TYPE    
df = df.with_columns(
    #cs.by_index(1).str.strip_prefix("chr").alias("Chr"),
    cs.by_index(1).alias("Chr"),
    cs.by_index(2).alias("Start"),
    cs.by_index(3).alias("End"),
    cs.by_index(4).alias("TYPE"),
    Strand = pl.lit(".")
    )

"""
df = df.with_columns(
        pl.when(
            #max and min should be above or equal 2
            (pl.col("Copy_Number").list.max() >= 2) &
            (pl.col("Copy_Number").list.min() >= 2)
        )
        .then(pl.lit("DUP"))
        .when(
            (pl.col("Copy_Number").list.max() < 2) &
            (pl.col("Copy_Number").list.min() < 2)
        )
        .then(pl.lit("DEL"))
        .otherwise(pl.lit("MIX"))
        .alias("CN_Type")
     )"""


(df.select(["Chr", "Start", "End", "TYPE", "Strand"])
   .unique(keep="any")
   .sort(by=["Chr", "Start", "End"])
   .sink_csv(output_file, separator = "\t", include_header = False))