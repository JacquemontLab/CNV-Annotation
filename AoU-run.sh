#!/bin/bash
~/nextflow run main.nf  --cnvs "gs://fc-secure-8bde5181-67bc-4770-96db-3b707c3f187f/output_cnvcalling/cnv_annotation/testing/penncnv_quantisnp_cnv_merged.tsv" \
--genome_version "GRCh38" \
--cohort_tag "AllOfUs" \
-c conf/aou.config \
-with-report \
-output-dir ~/cnv_annotation/DB-Builder \
-resume

