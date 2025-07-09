#!/usr/bin/env bash


###########################
#
# DB-Builder Cloud Set-Up
#
# Usage : bash cloud_setup.sh
# Description : For installing all necessary dependancies on fresh cloud VM. To be run from the setup directory. By default the program will install resources
#               (the vep cache and gnomad files) into the resources folder. Otherwise the location can be specified uing the -r flag.
#
# Requirements: - Requires a .bashrc in the $HOME
#               - curl
#               - git 
#               - perl
###########################

set -euo pipefail
IFS=$'\n\t'

RESOURCE_DIR=$(realpath ../resources)

# Function for checking existing installs
command_exists() {
    command -v "$1" &> /dev/null
}


print_usage() {
    echo "Usage: $0 [-r <resource_dir>]"
    echo
    echo "Options:"
    echo "  -r   Path to the resource directory (default: ../resources)"
    echo "  -h   Show this help message"
}

while getopts ":r:h" opt; do
    case ${opt} in
        r)
            RESOURCE_DIR=$(realpath "$OPTARG")
            ;;
        h)
            print_usage
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            print_usage
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            print_usage
            exit 1
            ;;
    esac
done

echo "Using resource directory: $RESOURCE_DIR"


#duckdb install 
if command_exists duckdb; then
    echo "skipping duckdb install"
else
    echo "Installing DuckDb..."
    curl https://install.duckdb.org | sh
    echo "Done!"
fi


#latest nextflow install
echo "Downloading latest Nextflow"

#We need at least Java v17
REQUIRED_MAJOR=17

check_java_version() {
    if ! command -v java &>/dev/null; then
        echo "Java is not installed."
        return 1
    fi

    JAVA_VERSION=$(java -version 2>&1 | awk -F[\".] '/version/ {print $2}')
    
    if [[ "$JAVA_VERSION" -ge "$REQUIRED_MAJOR" ]]; then
        echo "Java version $JAVA_VERSION is installed (OK)"
        return 0
    else
        echo "Java version $JAVA_VERSION is installed, but version $REQUIRED_MAJOR+ is required"
        return 1
    fi
}

if ! check_java_version; then
    curl -s https://get.sdkman.io | bash
    source ~/.bashrc
    sdk install java 17.0.10-tem
fi

#get nextflow and add to home.
curl -s https://get.nextflow.io | bash
chmod +x ~/nextflow
echo "New nextflow executable added at $HOME"

#Get VEP 113 
if ! command_exists vep;then
    echo "Adding VEP 113 to $HOME"
    git clone https://github.com/Ensembl/ensembl-vep.git $HOME
    cd $HOME/ensembl-vep  
    git checkout release/113
    perl INSTALL.pl
    cd -

    cat << 'EOF' >> $HOME/.bashrc
>>>>>> DB_BUILDER INSTALL >>>>>>>>
export PATH="$HOME/ensembl-vep:$PATH"
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
EOF
fi


#File dependancies

echo "Getting VEP cache"
cd $RESOURCE_DIR
curl -O https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz
tar xzf homo_sapiens_vep_113_GRCh38.tar.gz
curl -O https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh37.tar.gz
tar xzf homo_sapiens_vep_113_GRCh37.tar.gz

echo "Downloading extra gnomad resources"

#V4 Hg38
curl "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/exome_cnv/gnomad.v4.1.cnv.all.vcf.gz"  -o "gnomad.v4.1.cnv.all.vcf.gz" &&
                                                                        gunzip -c  "gnomad.v4.1.cnv.all.vcf.gz" |
                                                                        bgzip > "gnomad.v4.1.cnv.all.vcf.bgz"   &&
                                                                        tabix -p vcf "gnomad.v4.1.cnv.all.vcf.bgz"

curl "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz"  -o  "gnomad.v4.1.sv.sites.vcf.gz" &&
                                                                        gunzip -c  "gnomad.v4.1.sv.sites.vcf.gz" |
                                                                        bgzip > "gnomad.v4.1.sv.sites.vcf.bgz"   &&
                                                                        tabix -p vcf "gnomad.v4.1.sv.sites.vcf.bgz"

#V2 Hg19
curl "https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz" -o "gnomad_v2.1_sv.sites.vcf.gz" &&
                                                                        gunzip -c "gnomad_v2.1_sv.sites.vcf.gz" |
                                                                        bgzip > "gnomad.v2.1.sv.sites.vcf.bgz"   &&
                                                                        tabix -p vcf "gnomad.v2.1.sv.sites.vcf.bgz"
#Getting genetic Constraint metrics (pLi and LOEUF)
#https://gnomad.broadinstitute.org/help/constraint

curl "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv" -o "gnomad.v4.1.constraint_metrics.tsv"

echo "DONE"
