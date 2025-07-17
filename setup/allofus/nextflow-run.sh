#!/bin/bash

# Script to annotate a CNV TSV file locally.

# First argument is a .tsv file to annotate.
# It must contain at least the following columns: SampleID, Chr, Start, End, Type.
# Any additional columns present will be preserved and passed through.


NXF_VER=25.04.2
export NXF_OFFLINE=true

cnv_input_file=$1

nextflow run main.nf --cnvs "$cnv_input_file" \
    --genome_version "GRCh38" \
    --cohort_tag "AllOfUs" \
    -c setup/allofus/aou.config \
    -with-report report.html \
    -resume