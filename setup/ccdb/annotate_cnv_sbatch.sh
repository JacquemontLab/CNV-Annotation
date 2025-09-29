#!/bin/bash
# ------------------------------------------------------------------------------
# Script: annotate_cnv_sbatch.sh
#
# Description:
# This script submits a SLURM job to annotate a CNV TSV file using a Nextflow pipeline.
# It preserves all original columns and requires the input file to include at least the following:
#   - SampleID
#   - Chr
#   - Start
#   - End
#   - Type
#
# Usage:
#   sbatch annotate_cnv_sbatch.sh -i <CNV_TSV_FILE> -g <GENOME_VERSION> -c <COHORT_TAG>
#
# Inputs:
#   -d <GIT_DIR>        Path to the root of the repository containing `main.nf` and configs.
#   -i <CNV_TSV_FILE>   Path to a TSV file containing CNVs. Must include columns:
#                        SampleID, Chr, Start, End, Type. Additional columns are preserved.
#   -g <GENOME_VERSION>  Genome build (e.g., GRCh37, GRCh38) to use for annotation.
#   -c <COHORT_TAG>      Identifier for the cohort (used in annotation and output naming).
#
# Example usage:
#   sbatch annotate_cnv_sbatch.sh -i /path/to/input_cnvs.tsv -g GRCh38 -c MyCohort_Name -d /path/to/CNV-Annotation
# ------------------------------------------------------------------------------


#SBATCH --job-name=cnvDB_Buider    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=3500MB                     # Job memory request
#SBATCH --time=01:30:00               # Time limit hrs:min:sec
#SBATCH --output=cnv_annotation_%j.log   # Standard output and error log
#SBATCH --account=rrg-jacquese        #group account


# Default GIT_DIR
git_dir="/home/flben/links/projects/rrg-jacquese/LAB_WORKSPACE/SOFTWARE/Git_pipeline/CNV-Annotation/"

# ---------------------------
# Argument parsing
# ---------------------------
while getopts "i:g:c:d:" opt; do
  case $opt in
    i) cnv_input_file="$OPTARG" ;;
    g) genome_version="$OPTARG" ;;
    c) cohort_tag="$OPTARG" ;;
    d) git_dir="$OPTARG" ;;  # overrides default if provided
    *) echo "Usage: $0 -i <CNV_TSV_FILE> -g <GENOME_VERSION> -c <COHORT_TAG> -d <GIT_DIR>"
       exit 1 ;;
  esac
done


# Check required arguments
if [ -z "$cnv_input_file" ] || [ -z "$genome_version" ] || [ -z "$cohort_tag" ] || [ -z "$git_dir" ]; then
  echo "Error: Missing required arguments."
  echo "Usage: $0 -i <CNV_TSV_FILE> -g <GENOME_VERSION> -c <COHORT_TAG> -d <GIT_DIR>"
  exit 1
fi

# Check that main.nf exists in git_dir
if [ ! -f "${git_dir}/main.nf" ]; then
  echo "Error: Nextflow script '${git_dir}/main.nf' does not exist."
  exit 1
fi

# Get script directory local
# SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


# Load modules
module load java
module load python/3.13.2
module load nextflow
export NXF_OFFLINE=true

# If SLURM_TMPDIR is not set, fallback to system temp dir
TMPDIR=${SLURM_TMPDIR:-/tmp}

# Create and activate virtual environment
virtualenv --no-download "$TMPDIR/env"
source "$TMPDIR/env/bin/activate"

pip install --no-index --upgrade pip
pip install --no-index polars

# Run Nextflow pipeline
nextflow run "${git_dir}/main.nf" --cnvs "$cnv_input_file" \
    --genome_version "$genome_version" \
    --cohort_tag "$cohort_tag" \
    -c "${git_dir}/setup/ccdb/ccdb.config" \
    -with-report report.html \
    -resume
