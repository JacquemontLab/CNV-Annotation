#!/bin/bash

# Script to annotate a CNV TSV file locally.

# First argument is a .tsv file to annotate.
# It must contain at least the following columns: SampleID, Chr, Start, End, Type.
# Any additional columns present will be preserved and passed through.


export NXF_OFFLINE=true

cnv_input_file=$1


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

nextflow run ${SCRIPT_DIR}/../../main.nf --cnvs "$cnv_input_file" \
    --genome_version "GRCh38" \
    --cohort_tag "AllOfUs" \
    -c ${SCRIPT_DIR}/allofus.config \
    -with-report report.html \
    -resume