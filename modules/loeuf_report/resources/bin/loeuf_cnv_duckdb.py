#!/usr/bin/env python3
"""
LOEUF vs CNV Frequency Plot

This script generates a plot comparing LOEUF scores to CNV frequencies
across genes, using DuckDB for efficient data processing.

Usage:
    python loeuf_cnv_duckdb.py \
        -l path/to/loeuf_file.tsv \
        -c path/to/cnv_file.tsv \
        [-w window_size] \
        [-f overlap_column] \
        [-t overlap_threshold] \
        [-o output_plot.png]

Arguments:
    -l, --loeuf       Path to LOEUF file (TSV or Parquet)
    -c, --cnv         Path to CNV file (TSV or Parquet)
    -w, --window      Number of genes per window (default: 1000)
    -f, --overlap_col Column in CNV file for overlap filtering (optional)
    -t, --threshold   Threshold for overlap column (default: 0.5)
    -o, --output      Output plot filename (default: loeuf_cnv_plot.png)

Output:
    A PNG plot showing mean CNV observations per 1,000 individuals
    versus mean LOEUF scores across gene windows.
"""

import duckdb
import pyarrow
import pandas as pd
import polars as pl
import matplotlib.pyplot as plt
import argparse
import os
import sys

# -----------------------------
# Parse command-line arguments
# -----------------------------
parser = argparse.ArgumentParser(description="LOEUF vs CNV frequency plot using DuckDB")
parser.add_argument("-l", "--loeuf", required=True, help="Path to LOEUF file (TSV)")
parser.add_argument("-c", "--cnv", required=True, help="Path to CNV file (TSV)")
parser.add_argument("-w", "--window", type=int, default=1000, help="Window size [default 1000]")
parser.add_argument("-f", "--overlap_col", default=None, help="CNV overlap column (optional)")
parser.add_argument("-t", "--threshold", type=float, default=0.5, help="CNV overlap threshold [default 0.5]")
parser.add_argument("-o", "--output", default="loeuf_cnv_plot.png", help="Output plot file [default loeuf_cnv_plot.png]")
args = parser.parse_args()

# -----------------------------
# Check files exist
# -----------------------------
if not os.path.exists(args.loeuf):
    sys.exit(f"LOEUF file not found: {args.loeuf}")
if not os.path.exists(args.cnv):
    sys.exit(f"CNV file not found: {args.cnv}")

# -----------------------------
# Connect to DuckDB
# -----------------------------
con = duckdb.connect()

# -----------------------------
# Load LOEUF and CNV directly (Parquet or TSV)
# -----------------------------
def load_table(path):
    if path.endswith(".parquet"):
        return pl.scan_parquet(f"{path}")
    else:
        # default TSV
        return pl.scan_csv(f"{path}",separator = "\t", infer_schema_length=None)

loeuf = load_table(args.loeuf)
cnv_df = load_table(args.cnv)

# -----------------------------
# Check required columns
# -----------------------------
required_loeuf_cols = {"mane_select", "gene_id", "lof.oe_ci.upper"}
missing = required_loeuf_cols - set(loeuf.collect_schema().names())

if missing:
     sys.exit(f"LOEUF file missing columns: {', '.join(missing)}")

required_cnv_cols = {"SampleID", "Gene_ID"}
missing = required_cnv_cols - set(cnv_df.collect_schema().names())
if missing:
     sys.exit(f"CNV file missing columns: {', '.join(missing)}")

# -----------------------------
# Filter LOEUF
# -----------------------------

loeuf = loeuf.filter((pl.col("canonical") == True) &
                     (pl.col("gene_id").str.starts_with("ENS"))
                    )
#convert LOEUF column into float
loeuf = loeuf.with_columns(pl.col("lof.oe_ci.upper").replace("NA", None).cast(pl.Float64).alias("lof.oe_ci.upper"))

# Keep only rows where Exon_overlap > 0
cnv_df = cnv_df.filter(pl.col("Exon_Overlap") > 0).unique(subset=["SampleID", "Gene_ID"])

#pull number of unique Samples
nb_sample = (
    cnv_df
    .select(pl.col("SampleID").n_unique().alias("nb_sample"))
    .collect()
    .item()
)

# Register for SQL
con.register("loeuf", loeuf)
con.register("cnv_df", cnv_df)

# -----------------------------
# Function to compute stats in DuckDB
# -----------------------------
def compute_window_stats(filter_condition="1=1", group_name="All CNVs", genes_per_window=1000):
    query = f"""
    WITH gene_counts AS (
        SELECT Gene_ID, COUNT(*) AS freq
        FROM cnv_df
        WHERE {filter_condition}
          AND Gene_ID IN (SELECT gene_id FROM loeuf)
        GROUP BY Gene_ID
    ),
    merged AS (
        SELECT l.*, COALESCE(freq, 0) AS freq
        FROM loeuf l
        LEFT JOIN gene_counts g ON l.gene_id = g.Gene_ID
    ),
    ranked AS (
        SELECT *,
               ROW_NUMBER() OVER (ORDER BY CAST("lof.oe_ci.upper" AS DOUBLE)) AS rn
        FROM merged
    ),
    windowed AS (
        SELECT *,
           CAST(((rn - 1) / {genes_per_window}) + 1 AS INTEGER) AS window_id
        FROM ranked
    )
    SELECT window_id,
           AVG(CAST("lof.oe_ci.upper" AS DOUBLE)) AS mean_loeuf,
           AVG(freq / {nb_sample} * 1000) AS mean_freq,
           STDDEV(freq / {nb_sample} * 1000) / SQRT(COUNT(*)) AS sd_freq,
           COUNT(*) AS n_genes,
           SUM(CASE WHEN freq = 0 THEN 1 ELSE 0 END) AS n_zero_freq,
           '{group_name}' AS group_name
    FROM windowed
    GROUP BY window_id
    ORDER BY window_id
    """
    return con.execute(query).df()


# -----------------------------
# Compute stats
# -----------------------------

plot_data = compute_window_stats(genes_per_window=args.window)

# Stats for filtered CNVs (overlap threshold + problematic regions)
if args.overlap_col in cnv_df.columns:
    cond = f'"{args.overlap_col}" >= {args.threshold} AND "problematic_regions_Overlap" < 0.5'
    group_name = f'{args.overlap_col} >= {args.threshold} \n problematic_regions_Overlap < 0.5'
    filtered_stats = compute_window_stats(filter_condition=cond, group_name=group_name, genes_per_window=args.window)
    plot_data = pd.concat([plot_data, filtered_stats], ignore_index=True)


print(plot_data)

# -----------------------------
# Plot two figures in the same PNG
# -----------------------------
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(6, 10))  # two rows

# --- Top plot: all CNVs ---
ax = axes[0]
for group, df in plot_data.groupby("group_name"):
    ax.errorbar(df["mean_loeuf"], df["mean_freq"], yerr=df["sd_freq"],
                label=group, capsize=3, marker='o', linestyle='-')

ax.set_title(f"Window size {args.window}")
ax.set_xlabel("LOEUF")
ax.set_ylabel("Mean Observation per gene per 1k ind")
ax.set_xlim(0, 2)
ax.set_ylim(0, None)
ax.legend()


# --- Bottom plot: filtered CNVs ---
if args.overlap_col in cnv_df.columns:
    ax = axes[1]
    filtered_data = plot_data[plot_data["group_name"] == group_name]
    
    ax.errorbar(filtered_data["mean_loeuf"], filtered_data["mean_freq"], yerr=filtered_data["sd_freq"],
                label=group_name, capsize=3, marker='o', linestyle='-', color='tab:orange')
    ax.set_title(f"Window size {args.window}")
    ax.set_xlabel("LOEUF")
    ax.set_ylabel("Mean Obs per gene per 1k ind")
    ax.set_xlim(0, 2)
    ax.set_ylim(0, None)
    ax.legend()


plt.tight_layout()
plt.savefig(args.output, dpi=100)
print(f"Combined plot saved to: {args.output}")