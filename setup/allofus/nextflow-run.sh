#!/bin/bash
# ------------------------------------------------------------------------------
# Script: annotate_cnv_sbatch.sh
#
# Description:
#   The input TSV must contain at least:
#     - SampleID
#     - Chr
#     - Start
#     - End
#     - Type
#   Any additional columns are preserved.
#
# Usage:
#   annotate_cnv_sbatch.sh \
#     -i <CNV_TSV_FILE> \
#     -g <GENOME_VERSION> \
#     -c <DATASET_NAME> \
#     -d <GIT_DIR>
#
# Example:
#   sbatch annotate_cnv_sbatch.sh \
#     -i input_cnvs.tsv \
#     -g GRCh38 \
#     -c AllOfUs_tierv8_Array_GRCh38_Annotation \
#     -d /path/to/CNV-Annotation
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------------------------
while getopts "i:g:c:d:" opt; do
  case "${opt}" in
    i) cnv_input_file="${OPTARG}" ;;
    g) genome_version="${OPTARG}" ;;
    c) dataset_name="${OPTARG}" ;;
    d) git_dir="${OPTARG}" ;;
    *)
      echo "Usage: $0 -i <CNV_TSV_FILE> -g <GENOME_VERSION> -c <DATASET_NAME> -d <GIT_DIR>"
      exit 1
      ;;
  esac
done

# ------------------------------------------------------------------------------
# Validation
# ------------------------------------------------------------------------------
if [[ -z "${cnv_input_file}" || -z "${genome_version}" || -z "${dataset_name}" ]]; then
  echo "❌ Error: Missing required arguments."
  echo "Usage: $0 -i <CNV_TSV_FILE> -g <GENOME_VERSION> -c <DATASET_NAME> -d <GIT_DIR>"
  exit 1
fi

if [[ ! -f "${cnv_input_file}" ]]; then
  echo "❌ Error: CNV input file not found: ${cnv_input_file}"
  exit 1
fi

if [[ ! -f "${git_dir}/main.nf" ]]; then
  echo "❌ Error: main.nf not found in ${git_dir}"
  exit 1
fi

export NXF_OFFLINE=false

export CONDA_PKGS_DIRS=/home/jupyter/.conda/pkgs

nextflow run \
  "${git_dir}/main.nf" \
  --cnvs "${cnv_input_file}" \
  --genome_version "${genome_version}" \
  --dataset_name "${dataset_name}" \
  -c "${git_dir}/setup/allofus/allofus.config" \
  -with-report report.html \
  -resume
