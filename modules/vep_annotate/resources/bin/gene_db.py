#!/usr/bin/env python3
import polars as pl
import sys


"""
===============================================================================
Script Name   : gene_db.py
Author        : Benjamin Clark
Created       : 2025-04-30
Last Modified : 2025-04-30
Version       : 1.0.0
Python Version: 3.x
Description   : Reformats VEP CNV annotation for a CNV-GENE database.

Usage:
    python3 gene_db.py <in_file.parquet> <out_file.parquet>

Dependencies:
    - polars

===============================================================================
"""

__author__ = "Benjamin Clark"
__version__ = "1.0.1"



def main():
    """
    Script entry point. This function defines the execution order of the following functions.
    """
    # Lazy df creation
    df = pl.scan_parquet(sys.argv[1])

    # Initial cleaning that shouldn't be done in parallel, (yet?)
    df = make_null(df)
    df = make_exon_overlap(df)

    # Parallel processing 
    out = (
            df.pipe(make_transcript_overlap)
            .pipe(make_max_gnomad) #expects NULL values
            .pipe(make_CNV_ID)
            .pipe(make_canon_bool)
            .pipe(make_consequence_list)
          )

    # Outfile streaming to second positional argument
    out.rename({"Feature": "Transcript_ID","Gene": "Gene_ID"}).sink_parquet(sys.argv[2], compression="lz4")
    



def make_exon_overlap(df):
    """
    Function for producing a column containing the percentage overlap of CNVs over a gene's exon.

    !!!
    Currently expects the OverlapPC column to be a string.
    !!!

    Parameters:
        df (pl.DataFrame): Input Polars DataFrame with at least the 'EXON' column,
                           and the 'OverlapPC' column.

    Returns:
        pl.DataFrame: A new DataFrame with the computed 'exon_overlap' column added.
    """
    # Getting range of overlapped exons
    df = df \
        .with_columns([
        pl.when(pl.col("EXON").is_not_null())
        .then(
            pl.col('EXON')
            .str.split(by="/")
            .list.get(0)
            .alias('exon_range')
        ).otherwise(None)
    ])

    # Getting total num of exons
    df = df \
        .with_columns([
        pl.when(pl.col("EXON").is_not_null())
        .then(
            pl.col('EXON')
            .str.split(by="/")
            .list.get(1)
            .cast(pl.Float32)
        ).otherwise(None).alias('exon_sum'),

        # splitting exon ranges, we convert to "0" to avoid index errors later
        pl.when(pl.col("exon_range").is_not_null())
        .then(
            pl.col('exon_range')
            .str.split(by="-")
        ).otherwise(["0"])
        .cast(pl.List(pl.Float32))
        .alias('exon_range_split')
    ])

    # Adding exon overlap
    df = df \
        .with_columns(
        # check if only one exon is found overlapping
        pl.when((pl.col("exon_range_split").list.len() == 1) &
                (pl.col("exon_range_split").list.get(0) != 0.0)
                )  # the resulting overlap % is 1/#_of_exons
        .then(
            (1.0 / pl.col("exon_sum"))
        )  # check if there is an overlap with multiple exons
        .when(
            pl.col("exon_range_split").list.len() == 2
        )  #
        .then(
            (pl.col("exon_range_split").list.max() - pl.col("exon_range_split").list.min() + 1) /
            pl.col("exon_sum")
        )
        .when(
            pl.col("OverlapPC").str.contains("100")
        )
        .then(
            pl.lit(1.0)
        )
        .otherwise(None)
        .alias("Exon_Overlap")
    ).drop("exon_range",
           "exon_range_split",
           "exon_sum")
    return df


def make_null(df):
    """
    Replaces placeholder '-' values with nulls across all columns in the DataFrame,
    except for the 'CANONICAL' column.

    Parameters:
        df (pl.DataFrame): Input Polars DataFrame.

    Returns:
        pl.DataFrame: A copy of the DataFrame with '-' values replaced by nulls.
    """
    df = df \
        .with_columns([
        pl.when(pl.col(col) == '-')
        .then(None)
        .otherwise(pl.col(col))
        .alias(col)
        for col in [c for c in df.collect_schema().names()] 
    ])
    return df


def make_transcript_overlap(df):
    """
    Converts the 'OverlapPC' column from a percentage string or value to a float 
    representing the proportion of transcript base pair overlap.
    
    Parameters:
        df (pl.DataFrame): Input Polars DataFrame containing an 'OverlapPC' column 
                           with overlap percentages.

    Returns:
        pl.DataFrame: A copy of the DataFrame with the 'transcript_bp_overlap' column added
                      and the 'OverlapPC' column removed.
    """
    df = df.with_columns(
        (pl.col("OverlapPC").cast(pl.Float32) / 100 ).alias("Transcript_BP_Overlap")
        ).drop("OverlapPC")
    return df

def make_canon_bool(df):
    df = df.with_columns(
        pl.when(pl.col("Feature").is_not_null())
        .then(
            pl.col("CANONICAL").fill_null("")  # Replace null with empty string
            .str.contains("YES")               # Will now return False or True
        )
        .otherwise(None)
        .cast(pl.Boolean)
        .alias("CANONICAL")
    )
    return df

def make_consequence_list(df):
    df = df.with_columns(
        pl.col("Consequence").str.split(",").alias("Consequence")
    )
    return df


def make_max_gnomad(df):
    """
    Processes gnomAD allele frequency columns in the input DataFrame to compute the 
    maximum observed frequency per variant across all relevant gnomAD populations.

    !!!!!
    This function expects that all "-" from VEP in gnomad columns are first transformed
    to null.
    !!!!!

    Parameters:
        df (pl.DataFrame): Input Polars DataFrame with columns named like 'gnomad_AF_*'
                           containing comma-separated allele frequencies.

    Returns:
        pl.DataFrame: A copy of the DataFrame with the 'gnomad_max_freq' column added
                      and all original 'gnomad_AF_*' columns removed.
    """
    # Getting maximum value per list element
    df = df \
        .with_columns([
            # Replace nulls with a Float 0
            pl.when(
                pl.col(col).is_null()
            ) 
            .then(pl.lit(0.0, dtype=pl.Float32))
            # Rows containing string is split, cast to float and the max is extracted
            .otherwise(
                 pl.col(col)
                .str.split(",")
                .list.eval(pl.element().cast(pl.Float32))
                .list.max()
            )
            .alias(col)
            for col in df.collect_schema() if col.startswith("gnomad")]
        )
    # Getting maximum across columns (rowwise)
    df = df \
        .with_columns(
            pl.max_horizontal(pl.col("^gnomad_.*$"))
            .alias("Gnomad_Max_AF")
        ).drop(pl.col("^gnomad_.*$"))
    return df


def make_CNV_ID(df):
    """
    Generates a standardized 'CNV_ID' column by combining genomic location and CNV type 
    information from the 'Location' and 'Allele' columns.

    Parameters:
        df (pl.DataFrame): Input Polars DataFrame with 'Location' and 'Allele' columns.

    Returns:
        pl.DataFrame: A copy of the DataFrame with a new, first-position 'CNV_ID' column.
    """
    df = df.with_columns(
    # Extract first three letters
        pl.col("Allele").str.slice(0,3).str.to_uppercase()
        .alias("Allele")
    )
    # Concat location and CNV type into one column
    df = df.with_columns(
        (pl.col("Location")  + "_" +  pl.col("Allele"))
        .alias("CNV_ID")
    )
    # Replace new coord delimiter
    df = df.with_columns(
        pl.col("CNV_ID").str.replace_all(r"[:|-]", "_")
        .alias("CNV_ID")
    )
    # Reorder to first position, remove old id
    cols = ["CNV_ID"] + [col for col in df.columns if col not in ["CNV_ID", "#Uploaded_variation"]]
 
    return df.select(cols)

if __name__ == "__main__":
    main()
