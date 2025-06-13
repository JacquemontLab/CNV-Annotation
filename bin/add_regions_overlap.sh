#!/bin/bash

###############################################################################
# Script Name: add_regions_overlap.sh
# Description: This script prepares region-specific BED files (telomere, centromere,
#              and segmental duplications) from a genome regions file and computes
#              the overlap of CNVs with these regions. This is useful for annotating
#              CNVs that intersect unreliable genomic regions.
#
# Usage:
#   ./compute_overlap_regions_test.sh <input_cnv_file> <regions_file> <genome_version> <output_file>
#
# Arguments:
#   <input_cnv_file>   : Input CNV table to annotate (TSV format)
#                         Must contain columns: SampleID, Chr, Start, End (tab-separated).
#   <regions_file>     : Genome regions file (TSV) with telomere/centromere/segmentaldup data
#   <genome_version>   : Genome version string (e.g., GRCh37, GRCh38)
#   <output_file>      : Path to output file with added region overlap annotations
#
# Author: Florian Bénitière
# Date: April 2025
###############################################################################

set -euo pipefail

# --------------------- Check input arguments --------------------- #
if [ "$#" -ne 4 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: $0 <input_cnv_file> <regions_file> <genome_version> <output_file>"
    echo ""
    echo "Arguments:"
    echo "  <input_cnv_file>   : Path to input CNV file (TSV format)"
    echo "  <regions_file>     : Path to region annotation file (TSV) with columns: Chr, Start, End, Region_Type, Genome_Version"
    echo "  <genome_version>   : Genome version string (e.g., GRCh37, GRCh38)"
    echo "  <output_file>      : Path to output file with CNV region overlap annotations"
    exit 1
fi

# Assign arguments
input_cnv_file="$1"
regions_file="$2"
genome_version="$3"
output_file="$4"



# Create a temporary BED file for segmental duplications
segmentaldup_db=$(mktemp --suffix=.bed)
awk -v genome="$genome_version" 'BEGIN {
    OFS="\t"
    print "Chr", "Start", "End"
}
NR > 1 && $4 == "segmentaldup" && $5 == genome {
    print $1, $2, $3
}' "$regions_file" > "$segmentaldup_db"

# Create a temporary BED file for par regions
par_db=$(mktemp --suffix=.bed)
awk -v genome="$genome_version" 'BEGIN {
    OFS="\t"
    print "Chr", "Start", "End"
}
NR > 1 && ($4 == "PAR1" || $4 == "PAR2" || $4 == "XTR") && $5 == genome {
    print $1, $2, $3
}' "$regions_file" > "$par_db"



# Format string to pass to overlap computation script
regions_to_overlap=segmentaldup:$segmentaldup_db,par:$par_db"



# Annotate input CNV file with overlap metrics
compute_regions_overlap_fraction.sh "$input_cnv_file" "$regions_to_overlap" "$output_file"


# Clean up temporary BED files
rm -f  "$segmentaldup_db" "$par_db"
