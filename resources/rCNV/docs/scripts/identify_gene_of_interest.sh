#!/bin/bash
# ============================================
# identify_gene_of_interest.sh
#
# This script processes CNV regions for both GRCh37 and GRCh38.
# It:
#   1. Converts CNV BED regions into dummy VCF variants
#   2. Annotates them with VEP (Variant Effect Predictor)
#   3. Filters annotated results with DuckDB
#   4. Removes transcripts overlapping problematic genomic regions
#   5. Outputs a cleaned list of protein-coding transcripts per CNV
#
# Dependencies:
#   - Apptainer with VEP container
#   - DuckDB
#   - CNV-Caller repository (for overlap computation scripts and region definitions)
#
# Author: Florian Bénitière
# ============================================

# Loop over genome versions
for genome_version in GRCh37 GRCh38; do
    echo ">>> Processing $genome_version"

    ################ STEP 1: CREATE DUMMY VCF ################
    # Convert CNV BED file into a fake VCF file so that VEP can process it.
    # Each CNV is represented as a "DEL" variant with coordinates from BED.
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "DEL", ".", ".", ".", "."}' "cnv_regions_${genome_version}.bed" \
        | sort -k1,1 -k2,2n > "dummy_${genome_version}.vcf"

    ################ STEP 2: RUN VEP ################
    # Annotate the dummy VCF using Ensembl VEP inside an Apptainer container.
    # Options:
    #   --offline   : use local cache instead of Ensembl API
    #   --fork 64   : parallelize across 64 threads
    #   --overlaps  : report overlap with known features
    #   --biotype   : add gene biotype info
    #   --canonical : keep only canonical transcripts
    apptainer exec /lustre09/project/6008022/LAB_WORKSPACE/SOFTWARE/VEP/vep.sif \
        vep --input_file "dummy_${genome_version}.vcf" \
            --dir_cache /lustre09/project/6008022/LAB_WORKSPACE/SOFTWARE/VEP/cache \
            --offline \
            --fork 64 \
            --species homo_sapiens \
            --assembly $genome_version \
            --tab \
            --overlaps \
            --biotype \
            --canonical \
            --output_file "dummy_annotated_${genome_version}.txt" \
            --force_overwrite

    ################ STEP 3: FILTER WITH DUCKDB ################
    # Parse VEP output and retain only:
    #   - Canonical transcripts
    #   - Protein-coding genes
    #   - Variants overlapping > 50% with the region (OverlapPC > 50)
    # Output is written as a BED-like TSV (chr, start, end, gene, type).
    duckdb -c "
    COPY (
        SELECT DISTINCT
        b.Gene || '_' || 
        CAST(SPLIT_PART(SPLIT_PART(b.Location, ':', 2), '-', 1) AS BIGINT) || '_' ||  
        CAST(SPLIT_PART(SPLIT_PART(b.Location, ':', 2), '-', 2) AS BIGINT) AS SampleID,
        'chr' || SPLIT_PART(b.Location, ':', 1) AS Chr,
            t.Start,
            t.Stop,
            'DUP' AS Type,
        FROM read_csv_auto('dummy_annotated_${genome_version}.txt', delim='\t', header=true, all_varchar=true) b
        LEFT JOIN read_parquet('/home/flben/links/projects/rrg-jacquese/flben/Git/CNV-Annotation/resources/Transcript_Metadata/transcriptDB_${genome_version}.parquet') t
        ON b.Feature = t.Transcript_ID
        WHERE b.CANONICAL = 'YES'
        AND b.BIOTYPE = 'protein_coding'
        AND TRY_CAST(b.OverlapPC AS DOUBLE) > 50
    ) TO 'random_bed_${genome_version}_joined.tsv' (DELIMITER '\t', HEADER);
    "

    ################ STEP 3b: CHECK MATCHES ################
    # Count transcripts that match the database vs missing.
    duckdb -c "
    -- Count how many have matches vs missing
    SELECT
        COUNT(*) AS total_transcripts,
        SUM(CASE WHEN t.Transcript_ID IS NOT NULL THEN 1 ELSE 0 END) AS matched,
        SUM(CASE WHEN t.Transcript_ID IS NULL THEN 1 ELSE 0 END) AS not_matched
    FROM read_csv_auto('dummy_annotated_${genome_version}.txt', delim='\t', header=true, all_varchar=true) b
    LEFT JOIN read_parquet('/home/flben/links/projects/rrg-jacquese/flben/Git/CNV-Annotation/resources/Transcript_Metadata/transcriptDB_${genome_version}.parquet') t
    ON b.Feature = t.Transcript_ID
    WHERE b.CANONICAL = 'YES'
    AND b.BIOTYPE = 'protein_coding'
    AND TRY_CAST(b.OverlapPC AS DOUBLE) > 50;
    "


    ################ STEP 4: REMOVE PROBLEMATIC TRANSCRIPTS ################
    # Load predefined problematic genomic regions and filter transcripts overlapping them by >=50%.
    regions_file=/home/flben/links/projects/rrg-jacquese/flben/Git/CNV-Caller/resources/Genome_Regions/Genome_Regions_data.tsv
    problematicregions_db=$(mktemp --suffix=.bed)

    # Extract problematic regions for this specific genome version
    awk -v genome="$genome_version" 'BEGIN {
        OFS="\t"; print "Chr", "Start", "End"
    }
    NR > 1 && $4 == "problematic_regions" && $5 == genome {
        print $1, $2, $3
    }' "$regions_file" > "$problematicregions_db"

    # Prepare region string for overlap script
    regions_to_overlap="ProblematicRegions:$problematicregions_db"

    # Run overlap computation
    /home/flben/links/projects/rrg-jacquese/flben/Git/CNV-Caller/modules/merge_dataset_CNV/resources/bin/compute_regions_overlap_fraction.sh \
                 "random_bed_${genome_version}_joined.tsv" "$regions_to_overlap" "overlap_problematic_${genome_version}.tsv"

    ################ STEP 5: FINAL FILTER ################
    # Keep only transcripts with <50% overlap with problematic regions.
    # Output: list of protein-coding transcripts per CNV after removing problematic ones.
    duckdb -c "
    COPY (
        SELECT 
            REPLACE(Chr, 'chr', '') AS Chr, 
            SPLIT_PART(SampleID, '_', 2) AS Start_rCNV,
            SPLIT_PART(SampleID, '_', 3) AS End_rCNV,
            SPLIT_PART(SampleID, '_', 1) AS Gene
        FROM read_csv_auto('overlap_problematic_${genome_version}.tsv', delim='\t', header=true, all_varchar=true)
        WHERE TRY_CAST(ProblematicRegions_Overlap AS DOUBLE) < 0.5
    ) TO 'dummy_annotated_filtered_genes_${genome_version}.tsv' (DELIMITER '\t', HEADER);
    "

    echo ">>> Finished $genome_version"
done
