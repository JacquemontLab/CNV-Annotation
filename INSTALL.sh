#!/bin/bash

# ####################################
# CNV-Annotation Set-Up
#
# Usage: bash INSTALL.sh
#
# Description : For installing all necessary dependancies. By default the program will install resources
#               (the vep cache and gnomad files) into the resources folder. Otherwise the location can be specified uing the -r flag.
#
# Options:
#   -r <path> : Path to the resource directory (default: ./resources)
#
# Requirements:
#   - curl
#   - git
#   - perl
#   - .bashrc present in $HOME
# ####################################


set -e  # Exit immediately if a command exits with a non-zero status

#--- Logging function ---
log_step() {
    # Choose a symbol depending on the message
    case "$1" in
        STEP*) icon="🔹" ;;   # blue diamond for pipeline steps
        ERROR*) icon="❌" ;;  # red cross for errors
        WARN*) icon="⚠️ " ;;  # warning sign
        DONE*) icon="✅" ;;   # check mark for done
        *) icon="ℹ️ " ;;      # info icon
    esac
    echo -e "\n[$(date '+%Y-%m-%d %H:%M:%S')] $icon $1"
}

# --- Default Parameters ---
RESOURCE_DIR="$(git rev-parse --show-toplevel)/resources"

# --- Git-Variables
GIT_PROJECT="CNV-Annotation"


# --- Command-Line Argument Parser ---
print_usage() {
    echo "Usage: $0 [-r <resource_dir>]"
    echo "  -r   Path to the resource directory (default: ./resources)"
    echo "  -g   Genome Version for VEP Cache. Either GRCh38 or GRCh37"
    echo "  -h   Show this help message"
}

while getopts ":r:h:g:" opt; do
    case "${opt}" in
        r) RESOURCE_DIR="$(realpath "$OPTARG")" ;;
        g) GENOME_ASSEMBLY="$OPTARG" ;;
        h) print_usage && exit 0 ;;
        \?) echo "Invalid option: -$OPTARG" >&2; print_usage; exit 1 ;;
        :)  echo "Option -$OPTARG requires an argument." >&2; print_usage; exit 1 ;;
    esac
done

log_step "STEP Checking launch context"

# Validate GENOME_ASSEMBLY
if [[ "$GENOME_ASSEMBLY" != "GRCh38" && "$GENOME_ASSEMBLY" != "GRCh37" ]]; then
    log_step "ERROR Invalid genome version: $GENOME_ASSEMBLY (expected GRCh37 or GRCh38)"
    exit 1
fi

# ---  Checking if we are in the git repo if default resource path is being used ---
#check if default is being used
if [[ "$RESOURCE_DIR" == "$(git rev-parse --show-toplevel)/resources" ]]; then
        # check if we are in a git repo
        GIT_REPO=$(git rev-parse --show-toplevel 2>/dev/null)
        if [[ -n "$GIT_REPO" ]]; then
                repo_name=$(basename "$GIT_REPO")
        else
            log_step "ERROR Invalid launch context: git repository not found"
            exit 1
        fi
        # check if we are in the expected repo
        if [[ "$GIT_PROJECT" == "$repo_name" ]]; then
            log_step "DONE Default resource directory found at: $RESOURCE_DIR"
        else
            log_step "ERROR Invalid launch context: please launch script within the CNV-Annotation repository or provide -r"
            exit 1
        fi
fi

log_step "STEP Using resource directory: $RESOURCE_DIR"
mkdir -p "$RESOURCE_DIR"

# --- Function Definitions ---
add_to_path_once() {
    local line="$1"
    grep -qxF "$line" "$HOME/.bashrc" || echo "$line" >> "$HOME/.bashrc"
}


# --- Install Nextflow 25.10.2 if not present ---
log_step "STEP Checking Nextflow installation"
NF_REQUIRED_VERSION="25.10.2"

if command -v nextflow &> /dev/null; then
    NF_CURRENT_VERSION=$(nextflow -version | head -n 3 | grep version | awk '{print $2}')
    if [[ "$NF_CURRENT_VERSION" == "$NF_REQUIRED_VERSION" ]]; then
        log_step "DONE Nextflow $NF_REQUIRED_VERSION is already installed"
        INSTALL_NEXTFLOW=false
    else
        log_step "WARN Nextflow $NF_CURRENT_VERSION found, required $NF_REQUIRED_VERSION"
        INSTALL_NEXTFLOW=true
    fi
else
    log_step "WARN Nextflow not found"
    INSTALL_NEXTFLOW=true
fi

if [[ "${INSTALL_NEXTFLOW:-false}" == true ]]; then
    log_step "STEP Installing Nextflow $NF_REQUIRED_VERSION"
    curl -L -o "nextflow" https://github.com/nextflow-io/nextflow/releases/download/v25.10.2/nextflow
    #curl -s https://get.nextflow.io | bash
    mkdir -p "$HOME/bin"
    mv nextflow "$HOME/bin/"
    chmod +x "$HOME/bin/nextflow"   
    add_to_path_once 'export PATH="$HOME/bin:$PATH"'
    export PATH="$HOME/bin:$PATH"
    log_step "DONE Nextflow installed successfully"
fi


# --- Download VEP Cache ---
log_step "STEP Preparing VEP cache for $GENOME_ASSEMBLY"
pushd "$RESOURCE_DIR" > /dev/null

cache_file="homo_sapiens_vep_113_${GENOME_ASSEMBLY}.tar.gz"
cache_path="$RESOURCE_DIR/$cache_file"
cache_dir="$RESOURCE_DIR/homo_sapiens/113_${GENOME_ASSEMBLY}"

if [[ ! -d "$cache_dir" ]]; then
    log_step "STEP Downloading VEP cache ($GENOME_ASSEMBLY)"
    curl -o "$cache_path" "https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/${cache_file}"
    tar -xzf "$cache_path" -C "$RESOURCE_DIR"
    log_step "DONE VEP cache extracted"
else
    log_step "DONE VEP cache already exists at $cache_dir"
fi


# --- Download GnomAD Resources ---
log_step "STEP Downloading GnomAD CNV/SV resources"

# Create and enter ressources_gnomAD directory
mkdir -p ressources_gnomAD
pushd ressources_gnomAD > /dev/null

download_and_index_vcf() {
    local url="$1"
    local base=$(basename "$url")
    local bgz="${base/.vcf.gz/.vcf.bgz}"

    if [[ ! -f "$bgz" ]]; then
        log_step "STEP Downloading $(basename "$url")"
        curl -O "$url"
        gunzip -c "$base" | bgzip > "$bgz"
        tabix -p vcf "$bgz"
        log_step "DONE Indexed $bgz"
    else
        log_step "DONE $bgz already exists"
    fi
}

# GnomAD v4.1 (GRCh38)
download_and_index_vcf "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/exome_cnv/gnomad.v4.1.cnv.all.vcf.gz"
download_and_index_vcf "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz"

# GnomAD v2.1 (GRCh37)
download_and_index_vcf "https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz"

popd > /dev/null  # Exit ressources_gnomAD


# --- Constraint Metrics ---
log_step "STEP Downloading constraint metrics"

# Create and enter ressources_LOEUF directory
mkdir -p ressources_LOEUF
pushd ressources_LOEUF > /dev/null

if [[ ! -f "gnomad.v4.1.constraint_metrics.tsv" ]]; then
    curl -O "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
    log_step "DONE Constraint metrics downloaded"
else
    log_step "DONE Constraint metrics already present"
fi

popd > /dev/null  # Exit ressources_LOEUF
popd > /dev/null  # Exit RESOURCE_DIR

log_step "DONE All downloads complete"
log_step "DONE Setup finished successfully 🎉"