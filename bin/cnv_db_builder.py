#!/usr/bin/env python
import polars as pl
import sys

input = sys.argv[1]
output = sys.argv[2]


#defining explicit input datatypes, here I'm keeping a min of int32 since space
#should never be an issue
input_schema = pl.Schema(
    {
        "SampleID": pl.String,
        "Chr" : pl.String,
        "Start": pl.Int64,
        "End": pl.Int64,
        "Length": pl.Int32,
        "Copy_Number": pl.String,
        "Confidence_max": pl.Float32,
        "Num_Probes": pl.Int32,
        "Num_Merged_CNVs": pl.Int32,
        "QuantiSNP_Overlap": pl.Float32,
        "PennCNV_Overlap": pl.Float32,
        "Two_Algorithm_Overlap": pl.Float32,
        "telomere_Overlap": pl.Float32,
        "centromere_Overlap": pl.Float32,
        "segmentaldup_Overlap": pl.Float32
    }
)

#scan csv and cast to schema
df  = pl.scan_csv(input,
                  separator="\t",
                  schema=input_schema)
#Create Copy number column which is a list of variable small integers
df = df.with_columns(
    (pl.col("Copy_Number")
        .str.split(",")
        .list.eval(pl.element().cast(pl.Int32))
     ).alias("Copy_Number"))

#adding DEL,DUP,MIX flags
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
#Creating CNV_ID Column
df = df.with_columns(
        pl.concat_str([pl.col("Chr"),
                      pl.col("Start"),
                      pl.col("End"),
                      pl.col("CN_Type")],
                      separator = "_").alias("CNV_ID")
)

#Defining output order, putting IDs in front
order = (["CNV_ID", "SampleID"] +
         [col for col in df.collect_schema().names() if col not in ["CNV_ID", "SampleID"]])

df = df.select(order)
df.sink_parquet(output)
