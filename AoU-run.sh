#!/bin/bash
~/nextflow run AoU_main.nf  --cnvs "/home/jupyter/cnv_annotation/inputs/DB-Builder/test_v1_cnv_merged.tsv" \
--genome_version "GRCh38" \
--cohort_tag "AllOfUs" \
-c conf/aou.config \
-with-report \
-output-dir ~/cnv_annotation/DB-Builder

