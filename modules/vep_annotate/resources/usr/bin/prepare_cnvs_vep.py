#!/usr/bin/env python

import polars as pl
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]


df  = pl.scan_csv(input_file, separator = "\t").select(["Chr","Start","End", "Copy_Number"])

df = df.with_columns(
    (pl.col("Copy_Number")
        .str.split(",")
        .list.eval(pl.element().cast(pl.Int8))
     ).alias("Copy_Number"),

    (pl.col("Chr").str.strip_prefix("chr")).alias("Chr"),

    Strand = pl.lit(".")
    )

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
    )


(df.select(["Chr","Start","End", "CN_Type", "Strand"])
   .unique(keep="any")
   .sort(by=["Chr", "Start", "End"])
   .sink_csv(output_file, separator = "\t", include_header = False))