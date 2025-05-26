#!/bin/bash

###############################################################################
# Script Name: compute_cnv_region_overlap_fraction.sh
# Description: Calculates the fraction of overlap between CNV intervals and
#              multiple sets of genomic regions (BED files). Outputs a TSV
#              file with overlap fractions for each region set.
#
# Usage:
#   ./annotate_cnv_overlap_by_regions.sh <cnv_file[.gz]> <region_list> <output_file[.gz]>
#
#   <cnv_file[.gz]>     : CNV input file, gzipped or plain text.
#                         Must contain columns: SampleID, Chr, Start, End (tab-separated).
#   <region_list>       : Comma-separated list of region definitions in the form:
#                         "name1:file1.bed[.gz],name2:file2.bed[.gz],..."
#                         Each BED file must have: Chr, Start, End.
#   <output_file[.gz]>  : Output TSV file (gzipped if ending in .gz).
#
# Author: Florian Bénitière
# Date: April 2025
###############################################################################

# Check input arguments
if [ "$#" -ne 3 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: $0 <input_file[.gz]> <regions_to_overlap[.gz]> <output_file[.gz]>"
    echo ""
    echo "Arguments:"
    echo "  <input_file[.gz]>       : The input CNV file (either gzipped or uncompressed)."
    echo "                            Must contain at least the following columns:"
    echo "                            SampleID\tChr\tStart\tEnd"
    echo "  <regions_to_overlap> : A comma-separated list of regions with corresponding BED files."
    echo "                            Each entry should be in the format 'first_region_name:bed_file,second_region_name:bed_file'."
    echo "                            Each BED file must contain at least the following columns:"
    echo "                            Chr\tStart\tEnd"
    echo "  <output_file[.gz]>      : The TSV output file where the results will be saved (either gzipped or uncompressed)."
    exit 1
fi


# Assign arguments
input_file="$1"
regions_to_overlap="$2"
output_file="$3"


# Choose the appropriate decompression command
if [[ "$input_file" == *.gz ]]; then
    read_cmd="zcat"
    if ! command -v zcat &>/dev/null; then
        read_cmd="gunzip -c"
    fi
else
    read_cmd="cat"
fi

# Create temporary files
cnv_bed=$(mktemp)
cnv_annotated=$(mktemp)

# Save header line
header_update=$($read_cmd "$input_file" | head -n 1)

# --- Format CNV into BED format (Chr, Start, End, SampleID) ---s
$read_cmd "$input_file" | tail -n +2 | awk 'BEGIN{OFS="\t"}{ print $2, $3, $4, $1}' | sort -k1,1 -k2,2n > "$cnv_bed"

# --- Save original CNV data (SampleID,Chr,Start,End...) but format for later pasting ---
$read_cmd "$input_file" | tail -n +2 | awk 'BEGIN{OFS="\t"} {$1=$1","$2; $2=""; print}' - | sort -k1,1 -k2,2n > "$cnv_annotated"


# --- Process each region set (name:file.bed) ---
for item in ${regions_to_overlap//,/ }; do
    
    # Extract region name (before the colon)
    region=${item%%:*}
    # Extract BED file name (after the colon)
    tsvfile=${item##*:}

    # Determine decompression method for region file
    if [[ "$tsvfile" == *.gz ]]; then
        read_cmd="zcat"
        if ! command -v zcat &>/dev/null; then
            read_cmd="gunzip -c"
        fi
    else
        read_cmd="cat"
    fi

    # --- Sort BED file ---
    region_bed_sorted=$(mktemp)
    $read_cmd "$tsvfile" | tail -n +2 | sort -k1,1 -k2,2n > "$region_bed_sorted"
    
    # --- Merge overlapping regions ---
    merged_region_bed=$(mktemp)
    bedtools merge -i "$region_bed_sorted" > "$merged_region_bed"

    # --- Intersect CNVs with merged regions, get overlap sizes ---
    overlap_bed=$(mktemp)
    bedtools intersect -a "$cnv_bed" -b "$merged_region_bed" -wao > "$overlap_bed"

    # --- Sum overlaps per CNV, calculate overlap fraction ---
    overlap_summary=$(mktemp)
    awk 'BEGIN{OFS="\t"} {$1=$4","$1; $4=""; print}' "$overlap_bed" | sort -k1,1 -k2,2n | bedtools merge -i - -c 8 -o sum > "$overlap_summary"

    # Merge overlapping intervals and compute fraction
    temp_frac_overlap_bed=$(mktemp)
    awk -F'\t' 'BEGIN{OFS="\t"} {
            size = $3 - $2;
            overlap = $NF;
            frac = overlap / size ;
            print $0, frac
        }' "$overlap_summary" > "$temp_frac_overlap_bed"

    # --- Append the fraction to the original CNV data ---
    updated_cnv_annotated=$(mktemp)
    paste "$cnv_annotated" <(cut -f5 "$temp_frac_overlap_bed") > "$updated_cnv_annotated"
    mv "$updated_cnv_annotated" "$cnv_annotated"

    # Update header
    header_update="$header_update\t${region}_Overlap"

    # Clean up per-region temp files
    rm -f "$region_bed_sorted" "$merged_region_bed" "$overlap_bed" "$overlap_summary"
done

# Update the final output with the new header and the results
# - Add the region names to the header
# - Merge the CNV data with the newly calculated overlap fractions
final_output=$(mktemp)
( echo -e "$header_update" && \
awk -F'\t' 'BEGIN{OFS="\t"} {split($1, arr, ","); $1=arr[1]; $2=arr[2]; print $0}' "$cnv_annotated" ) | {
    if [[ "$output_file" == *.gz ]]; then
        gzip > "$final_output"
    else
        cat > "$final_output"
    fi
}

mv "$final_output" "$output_file"

# --- Clean up general temporary files ---
rm -f "$cnv_bed" "$cnv_annotated"