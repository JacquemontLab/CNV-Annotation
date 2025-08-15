#!/bin/bash

# Script to annotate a CNV TSV file in a job.

# First argument is a .tsv file to annotate.
# It must contain at least the following columns: SampleID, Chr, Start, End, Type.
# Any additional columns present will be preserved and passed through.

#SBATCH --job-name=cnvDB_Buider    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=3500MB                     # Job memory request
#SBATCH --time=01:00:00               # Time limit hrs:min:sec
#SBATCH --output=cnvDB_Buider_%j.log   # Standard output and error log
#SBATCH --account=rrg-jacquese        #group account


module load java
module load python/3.13.2
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index polars


module load nextflow
export NXF_OFFLINE=true

cnv_input_file=$1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

nextflow run ${SCRIPT_DIR}/../../main.nf --cnvs "$cnv_input_file" \
    --genome_version "GRCh37" \
    --cohort_tag "cnvDB_Buider_GRCh37" \
    -c ${SCRIPT_DIR}/ccdb.config \
    -with-report report.html \
    -resume
