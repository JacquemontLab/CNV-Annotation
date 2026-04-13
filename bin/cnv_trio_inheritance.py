#!/usr/bin/env python3
# =============================================================================
# CNV Inheritance Analysis Pipeline
# Author: Florian Bénitière and Benjamin Clark
# Date: 22/07/2025
# Implementation: Python (Polars)
#
# Description:
# This pipeline analyzes the inheritance of CNVs (Copy Number Variants) 
# in trio datasets. It filters CNVs based on pedigree information, 
# constructs parent–child mappings, performs overlap checks with bedtools, 
# and outputs an annotated child CNV table.
#
# Requirements:
# - Python packages: polars, pyarrow
# - External tools: zcat, awk, bedtools
#
# Input:
# - CNV file (.tsv.gz or .parquet file)
# - Pedigree file (tab-separated: SampleID, FatherID, MotherID)
#
# Output:
# - Annotated child CNV BED and TSV files
# =============================================================================


import polars as pl
import os
import argparse
from pathlib import Path
import subprocess
pl.enable_string_cache()

# ------------------- Parse Arguments -------------------

parser = argparse.ArgumentParser(
    description="CNV Inheritance Analysis Pipeline (Polars implementation)"
)
parser.add_argument(
    "-c", "--cnv", required=True,
    help="Path to merged CNV file (.tsv.gz OR .parquet) (SampleID, Chr, Start, End)"
)
parser.add_argument(
    "-p", "--pedigree", required=True,
    help="Path to pedigree file (.tsv) (SampleID, FatherID, MotherID)"
)
parser.add_argument(
    "-o", "--output", default="annotated_child_cnv.tsv",
    help="Path to the output child CNVs"
)
parser.add_argument(
    "-t", "--type_col", default=None,
    help="Column name of CNV type (DEL, DUP, MIX etc...) for per-type annotation; leave empty to skip"
)
parser.add_argument(
    "-f", "--overlap", default="0.5",
    help="Comma-separated list of reciprocal overlap fractions (e.g., 0.5,0.1)"
)

args = parser.parse_args()

FILE_CNV = args.cnv
PEDIGREE_FILE = args.pedigree
OUTPUT = args.output
TYPE_COL = args.type_col
OVERLAPS = [float(x) for x in args.overlap.split(",")]

# Detect file type (TSV vs Parquet)
tsv_mode = True
suffix = Path(FILE_CNV).suffixes
if '.parquet' in suffix:
    tsv_mode = False
elif ".tsv" not in suffix:
    print("Invalid input type or malformed extension")
    sys.exit(1)


# ------------------- Load Pedigree -------------------
family_info = pl.scan_csv(PEDIGREE_FILE, separator="\t",infer_schema_length=1000000)
family_info = family_info.with_columns([
    pl.col("SampleID").cast(pl.Utf8),
    pl.col("FatherID").cast(pl.Utf8),
    pl.col("MotherID").cast(pl.Utf8)
])
family_info = family_info.select(["SampleID", "FatherID", "MotherID"])

# ------------------- Load CNVs -------------------
if TYPE_COL:
    if tsv_mode:
        df_cnv = pl.scan_csv(FILE_CNV, 
                             separator = "\t", 
                             schema_overrides={f"{TYPE_COL}": pl.Categorical, 
                                           "Chr": pl.Categorical,
                                           "SampleID": pl.Utf8},infer_schema_length=1000000)
        print(f"Loading {FILE_CNV} in tsv mode and using {TYPE_COL} as type column")
    #casting to categorical for performance 
    else:
        df_cnv = pl.scan_parquet(FILE_CNV)
        df_cnv = df_cnv.with_columns(pl.col(f"{TYPE_COL}").cast(pl.Categorical),
                                     pl.col("Chr").cast(pl.Categorical),
                                     pl.col("SampleID").cast(pl.Utf8))
        print(f"Loading {FILE_CNV} in parquet mode and using {TYPE_COL} as type column")
else:
    if tsv_mode:
        print(f"Loading {FILE_CNV} in tsv mode without type column")
        df_cnv = pl.scan_csv(FILE_CNV, 
                             separator = "\t", 
                             schema_overrides={"Chr": pl.Categorical,
                                               "SampleID": pl.Utf8},infer_schema_length=1000000)
    else:
        print(f"Loading {FILE_CNV} in parquet mode and without type column")
        df_cnv = pl.scan_parquet(FILE_CNV)
        df_cnv = df_cnv.with_columns(pl.col("Chr").cast(pl.Categorical))
        df_cnv = df_cnv.with_columns(pl.col("SampleID").cast(pl.Utf8))

print("Filtering for trios that are completely found in the cnv table...")

# keep only trios that have CNVs for every member 
sample_ids = df_cnv.select(pl.col("SampleID")).unique().collect().to_series()
family_trios =  family_info.filter(
    pl.col("SampleID").is_in(sample_ids.implode()) & #adding implodes as per https://github.com/pola-rs/polars/pull/22178
    pl.col("FatherID").is_in(sample_ids.implode()) &
    pl.col("MotherID").is_in(sample_ids.implode())
)

print("Building bed files...")

# ------------------- Subset Child and Parent CNVs -------------------
# Child CNVs (present in pedigree)
df_child = df_cnv.join(family_trios, on = "SampleID", how = "inner")

# Parent CNVs
df_father = df_cnv.join(family_trios, left_on = "SampleID", 
                                     right_on = "FatherID", 
                                     suffix = ".child", 
                                     how = "inner")
#mother cnvs
df_mother = df_cnv.join(family_trios, left_on = "SampleID", 
                                     right_on = "MotherID", 
                                     suffix = ".child", 
                                     how = "inner")


# Reorder columns and drop redundant IDs
standard_col = ["SampleID.child","Chr","Start","End"] 
f_order = standard_col + [col for col in df_father.collect_schema().names() if col not in standard_col]
m_order = standard_col + [col for col in df_mother.collect_schema().names() if col not in standard_col]
df_parents = pl.concat([df_father.select(f_order).drop("MotherID","SampleID"), 
                        df_mother.select(m_order).drop("FatherID", "SampleID")])

# ------------------- Build BED IDs for bedtools -------------------
if TYPE_COL is not None:
    # Parents
    parents_out = df_parents.with_columns(pl.concat_str([pl.col("SampleID.child"),
                                                         pl.lit("_"),
                                                         pl.col("Chr"), 
                                                         pl.lit("_"), 
                                                         pl.col(f"{TYPE_COL}")])
                                          .alias("SampleID")).select("SampleID","Start","End")
                    
    # Children                      
    df_child = df_child.with_columns(
        (pl.col("SampleID") + "_" + pl.col("Chr").cast(pl.Utf8) + "_" + pl.col(TYPE_COL) + "_" + pl.col("Start").cast(pl.Utf8) + "_" + pl.col("End").cast(pl.Utf8) ).alias("ID")
    )
    
    child_out = df_child.with_columns(pl.concat_str([pl.col("SampleID"),
                                                     pl.lit("_"),
                                                     pl.col("Chr"), 
                                                     pl.lit("_"), 
                                                     pl.col(f"{TYPE_COL}")])
                                          .alias("SampleID")).select("SampleID","Start","End")

else:

    parents_out = df_parents.with_columns(pl.concat_str([pl.col("SampleID.child"),
                                                         pl.lit("_"),
                                                         pl.col("Chr"), 
                                                         pl.lit("_")
                                                        ])
                                          .alias("SampleID")).select("SampleID","Start","End")

    df_child = df_child.with_columns(
        (pl.col("SampleID") + "_" + pl.col("Chr").cast(pl.Utf8) + "_" + pl.col("Start").cast(pl.Utf8) + "_" + pl.col("End").cast(pl.Utf8) ).alias("ID")
    )

    child_out = df_child.with_columns(pl.concat_str([pl.col("SampleID"),
                                                     pl.lit("_"),
                                                     pl.col("Chr"), 
                                                     pl.lit("_") 
                                                     ])
                                      .alias("SampleID")).select("SampleID","Start","End")


# ------------------- Write BED files -------------------
child_bed = f"child_cnv.bed"
parent_bed = f"parents_cnv.bed"

child_out.sink_csv(child_bed, separator = "\t", include_header = False)
parents_out.sink_csv(parent_bed, separator = "\t", include_header = False)


# ------------------- Run bedtools intersect -------------------
for ovlap in OVERLAPS:
    bash_cmd = f"bedtools intersect -a {child_bed} -b {parent_bed} -f {ovlap} -r -wa -u> intersect_ovlap{ovlap}.bed;"
    print("Running command:\n", bash_cmd)
    subprocess.run(bash_cmd, shell=True, check=True)

print("Formatting intersections and combining results...")

# ------------------- Annotate Child CNVs with Overlaps -------------------
for ovlap in OVERLAPS:
        intersects = pl.scan_csv(f"intersect_ovlap{ovlap}.bed", separator = "\t", has_header = False,infer_schema_length=1000000)
        
        # Build ID column for matching
        intersects = intersects.with_columns(
         (pl.col("column_1") + "_" + pl.col("column_2").cast(pl.Utf8) + "_" + pl.col("column_3").cast(pl.Utf8) ).alias("ID")
          )
          
    # Join with child CNVs and mark True/False
        df_child = df_child.join(
              intersects.select("ID").with_columns(
                  pl.lit(True).alias(f"Observed_in_Parent_{ovlap}")
              ),
              on="ID",
              how="left"
          ).with_columns(
              pl.col(f"Observed_in_Parent_{ovlap}").fill_null(False)
          )


# Drop temporary ID column
df_child = df_child.drop("ID")

print("Pipeline completed successfully.")


# ------------------- Write Final Annotated CNVs -------------------
df_child.drop("FatherID", "MotherID").sink_csv(f"{OUTPUT}", separator = "\t")
print(f"💾 Results written to {OUTPUT}")


# ------------------- Cleanup Temporary Files -------------------
files_to_remove = [
    parent_bed,
    child_bed
] + [f"intersect_ovlap{ov}.bed" for ov in OVERLAPS] 

for f in files_to_remove:
    try:
        os.remove(f)
        print(f"🗑️ Removed {f}")
    except FileNotFoundError:
        print(f"⚠️ Skipped {f} (not found)")