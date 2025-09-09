#!/usr/bin/env bash
# -*- coding: utf-8 -*-

# ####################################
# DB-Builder Cloud Set-Up Script
#
# Usage: bash cloud_setup.sh
#
# Description:
# Description : For installing all necessary dependancies on fresh cloud VM. By default the program will install resources
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

IFS=$'\n\t'


# --- Default Parameters ---
RESOURCE_DIR="$(git rev-parse --show-toplevel)/resources"


# --- Git-Variables
GIT_PROJECT="CNV-DB-Builder"


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

# ---  Checking if we are in the git repo if default resource path is being used ---

#check if default is being used
if [[ "$RESOURCE_DIR" == "$(git rev-parse --show-toplevel)/resources" ]]; then
        #check if we are in a git repo
        GIT_REPO=$(git rev-parse --show-toplevel 2>/dev/null)
        if [[ -n "$GIT_REPO" ]]; then
                repo_name=$(basename "$GIT_REPO")
        else
                echo "Invalid launch context: git repository not found"
                exit 1
        fi
        #check if we are in the expected repo
        if [[ "$GIT_PROJECT" == "$repo_name" ]]; then
                echo "Default resource directory found at: $RESOURCE_DIR"
        else
                echo "Invalid launch context: please launch script within the CNV-DB-Builder repository or"
                echo " provide a path to a resource directory (-r /path/to/dir)"
                exit 1
        fi
fi


# Validate GENOME_ASSEMBLY
if [[ "$GENOME_ASSEMBLY" != "GRCh38" && "$GENOME_ASSEMBLY" != "GRCh37" ]]; then
    echo "Invalid genome version: $GENOME_ASSEMBLY"
    exit 1
fi



echo "📁 Using resource directory: $RESOURCE_DIR"
mkdir -p "$RESOURCE_DIR"

# --- Function Definitions ---
command_exists() {
    command -v "$1" &>/dev/null
}


check_java_version() {
    if ! command_exists java; then
        echo "❌ Java is not installed."
        return 1
    fi
    local version
    version=$(java -version 2>&1 | awk -F[\".] '/version/ {print $2}')
    if [[ "$version" -ge 17 ]]; then
        echo "✅ Java version $version detected"
        return 0
    else
        echo "⚠️ Java version $version is too old, need 17+"
        return 1
    fi
}

add_to_path_once() {
    local line="$1"
    grep -qxF "$line" "$HOME/.bashrc" || echo "$line" >> "$HOME/.bashrc"
}

# --- Install DuckDB ---
if command_exists duckdb; then
    echo "✅ DuckDB already installed."
else
    echo "🔧 Installing DuckDB..."
    curl -fsSL https://install.duckdb.org | sh
    add_to_path_once 'export PATH="$HOME/.duckdb/cli/latest/:$PATH"'
fi

# --- Install Java (if needed) ---
if ! check_java_version; then
    echo "🔧 Installing Java 17 via SDKMAN..."
    curl -s https://get.sdkman.io | bash
    source "$HOME/.bashrc"
    source "$HOME/.zshrc"
    sdk install java 17.0.10-tem

fi

# --- Install Python Packages ---
echo "⬇️  Upgrading pip and installing Python packages in virtualenv"
ENV_DIR="$(git rev-parse --show-toplevel)/env"
mkdir -p "$ENV_DIR"
python3 -m venv "$ENV_DIR/db-builder-env"
source "$ENV_DIR/db-builder-env/bin/activate"
pip install --upgrade pip
pip install duckdb pandas matplotlib tqdm polars pyarrow
echo "✅ Python virtual environment installed in $ENV_DIR/db-builder-env. Be sure to source before run."
deactivate

# --- Install Nextflow ---
if command_exists nextflow; then
    echo "Nextflow command detected. Skipping installation..."
else
    echo "⬇️  Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    mkdir -p "$HOME/bin"
    mv nextflow "$HOME/bin/"
    chmod +x "$HOME/bin/nextflow"   
    add_to_path_once 'export PATH="$HOME/bin:$PATH"'
    export PATH="$HOME/bin:$PATH"
    echo "✅ Nextflow installed: $($HOME/bin/nextflow -version)"
fi
# --- Install VEP 113 ---
# VEP_DIR="$HOME/ensembl-vep"
# if ! command_exists vep; then
#     echo "⬇️  Cloning VEP 113 to $VEP_DIR"
#     if [[ -d "$VEP_DIR" && -n "$(ls -A "$VEP_DIR")" ]]; then
#         echo "⚠️  $VEP_DIR exists and is not empty. Skipping VEP clone."
#     else
#         git clone https://github.com/Ensembl/ensembl-vep.git "$VEP_DIR"
#         pushd "$VEP_DIR" > /dev/null
#         git checkout release/113
#         yes "n" |  perl INSTALL.pl
#         popd > /dev/null
#         add_to_path_once 'export PATH="$HOME/ensembl-vep:$PATH"'
#         export PATH="$HOME/ensembl-vep:$PATH"
#         echo "✅ VEP 113 installed."
#     fi
# else
#     echo "✅ VEP already installed."
#fi

# --- Download VEP Cache ---
echo "⬇️  Downloading VEP cache files..."
pushd "$RESOURCE_DIR" > /dev/null


cache_file="homo_sapiens_vep_113_${GENOME_ASSEMBLY}.tar.gz"
cache_path="$RESOURCE_DIR/$cache_file"
cache_dir="$RESOURCE_DIR/homo_sapiens/113_${GENOME_ASSEMBLY}"

if [[ ! -d "$cache_dir" ]]; then
    echo "⬇️  Downloading VEP cache for $GENOME_ASSEMBLY..."
    curl -o "$cache_path" "https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/${cache_file}"
    tar -xzf "$cache_path" -C "$RESOURCE_DIR"
else
    echo "✅ VEP cache for${GENOME_ASSEMBLY} already exists at $cache_dir"
fi



# --- Download GnomAD Resources ---
echo "⬇️  Downloading GnomAD CNV/SV files..."

# Create and enter ressources_gnomAD directory
mkdir -p ressources_gnomAD
pushd ressources_gnomAD > /dev/null

download_and_index_vcf() {
    local url="$1"
    local base=$(basename "$url")
    local bgz="${base/.vcf.gz/.vcf.bgz}"

    if [[ ! -f "$bgz" ]]; then
        curl -O "$url"
        gunzip -c "$base" | bgzip > "$bgz"
        tabix -p vcf "$bgz"
    else
        echo "✅ $bgz already exists."
    fi
}


# GnomAD v4.1 (GRCh38)
download_and_index_vcf "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/exome_cnv/gnomad.v4.1.cnv.all.vcf.gz"
download_and_index_vcf "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz"

# GnomAD v2.1 (GRCh37)
download_and_index_vcf "https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz"

popd > /dev/null  # Exit ressources_gnomAD


# --- Constraint Metrics ---
echo "⬇️  Downloading Constraint Metrics..."

# Create and enter ressources_LOEUF directory
mkdir -p ressources_LOEUF
pushd ressources_LOEUF > /dev/null

if [[ ! -f "gnomad.v4.1.constraint_metrics.tsv" ]]; then
    curl -O "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
else
    echo "✅ Constraint metrics already downloaded."
fi

popd > /dev/null  # Exit ressources_LOEUF
popd > /dev/null  # Exit RESOURCE_DIR

echo "✅ All downloads complete."
echo "🎉 Setup finished successfully!"