#!/usr/bin/env nextflow 



// This process formats CNV input files for VEP (Variant Effect Predictor) annotation.
// It extracts unique CNV coordinates to reduce redundant queries and speed up VEP processing.
process identifyUniqCNV {
    label 'quick'
    
    input:
    path cnvs 

    output:
    path "uniq_cnvs.bed"

    script:
    """
    prepare_cnvs_vep.py ${cnvs} "uniq_cnvs.bed"
    """
}

// This process annotates CNVs by overlapping them with genomic regions,
// then formats the output to include a unique CNV identifier for downstream analysis.

process computeOverlapRegion {    
    label 'quick'

    input:
    path uniq_cnvs
    val genomic_regions
    path regions_file 

    output:
    path "CNVs_overlap_region_with_CNV_ID.tsv"

    script:
    """
    # Add a SampleID column with value 'uniq' to the CNVs input file required for the script to run
    awk -v sample="uniq" \
    'BEGIN{OFS="\\t"; print "SampleID\\tChr\\tStart\\tEnd\\tType\\tStrand"} {print sample, \$0}' ${uniq_cnvs} > uniq_cnvs_with_sample.tsv

    # Run custom script to add overlap information between CNVs and genomic regions
    add_regions_overlap.sh uniq_cnvs_with_sample.tsv ${regions_file} ${genomic_regions} "CNVs_with_genomic_regions_overlap.tsv"

    # Format overlap output to add a CNV_ID column for unique CNV identification and keep only _Overlap columns
    format_overlap.sh "CNVs_with_genomic_regions_overlap.tsv" "CNVs_overlap_region_with_CNV_ID.tsv"
    """
}


// Annotating CNVs for gene disruptions using VEP (Variant Effect Predictor)
// Requires Singularity/Apptainer container for reproducibility and environment consistency

// Description: Runs VEP on unique CNV coordinates (GRCh38 assembly)
process VEP_GRCh38 {

    label 'vep'
    input:
    path uniq_cnvs
    path vep_cache
    path gnomad_sv


    output:
    path "vep_out.tsv", emit : results
    path "*html", emit : summary
    path "vep_comments.txt", emit : comments
    

    script:
    """
    tabix -p vcf ${gnomad_sv}

    vep -i ${uniq_cnvs} -o vep_out.tsv\
    -cache\
    --tab\
    --dir_cache ${vep_cache}\
    --offline\
    --force_overwrite\
    --numbers\
    --fork ${task.cpus}\
    --overlaps\
    --canonical\
    --max_sv_size 20000000\
    --verbose\
    --mane\
    --assembly GRCh38 \
    --custom file="./${gnomad_sv}",short_name=gnomad,format=VCF,reciprocal=1,overlap_cutoff=70,same_type=1,fields=AF_nfe%AF_afr%AF_amr%AF_fin%AF_sas%AF_eas%AF_asj \
    --fields "Uploaded_variation,Location,Allele,Gene,Feature,Consequence,CANONICAL,MANE,EXON,INTRON,OverlapPC,gnomad_AF_nfe,gnomad_AF_afr,gnomad_AF_amr,gnomad_AF_fin,gnomad_AF_sas,gnomad_AF_eas,gnomad_AF_asj"


    grep -E '^\\s*#' vep_out.tsv > vep_comments.txt
    """
}

// Description: Runs VEP on unique CNV coordinates (GRCh37 assembly)
process VEP_GRCh37 {
    label 'vep'
    
    input:
    path uniq_cnvs
    path vep_cache
    path gnomad_sv


    output:
    path "vep_out.tsv", emit : results
    path "*html", emit : summary
    path "vep_comments.txt", emit : comments
    

    script:
    """
    tabix -p vcf ${gnomad_sv}

    vep -i ${uniq_cnvs} -o vep_out.tsv\
    -cache\
    --tab\
    --dir_cache ${vep_cache}\
    --offline\
    --force_overwrite\
    --numbers\
    --fork ${task.cpus}\
    --overlaps\
    --canonical\
    --max_sv_size 20000000\
    --verbose\
    --assembly GRCh37 \
    --custom file="./${gnomad_sv}",short_name=gnomad,format=VCF,reciprocal=1,overlap_cutoff=70,same_type=1,fields=AFR_AF%AMR_AF%EAS_AF%EUR_AF \
    --fields "Uploaded_variation,Location,Allele,Gene,Feature,Consequence,CANONICAL,MANE,EXON,INTRON,OverlapPC,gnomad_AFR_AF,gnomad_AMR_AF,gnomad_EAS_AF,gnomad_EUR_AF"


    grep -E '^\\s*#' vep_out.tsv > vep_comments.txt
    """
}


// Process to build a gene database by integrating VEP annotations, gnomAD constraints variable LOEUF, 
// and transcript metadata. Outputs a compressed Parquet file with genome version metadata.

process buildGeneDB {

    input:
    path vep_out
    path gnomad_constraints
    path transcript_metadata
    path region_overlap
    val genome_version

    output:
    path "geneDB.parquet"

    script:
    """
    # First transform large VEP output to parquet
    duckdb -c "COPY (SELECT * FROM read_csv(${vep_out}, delim = '\\t')) 
               TO 'tmp_db.parquet' (FORMAT 'PARQUET', CODEC 'ZSTD');"

    # Formatting output
    gene_db.py tmp_db.parquet tmp_formatted.parquet

    # Adding gnomad_constraints file via right join on geneDB using gene_IDs
    duckdb -c "COPY ( 
                        SELECT geneDB.*, CAST(NULLIF(\\"lof.oe_ci.upper\\", 'NA') AS DOUBLE) AS LOEUF
                        FROM read_csv(\\"${gnomad_constraints}\\", delim = '\\t') AS gnomad
                        RIGHT JOIN (SELECT * FROM read_parquet('tmp_formatted.parquet')) AS geneDB ON geneDB.Transcript_ID = gnomad.transcript
                        
        ) TO "tmp_gene_constraints.parquet" (FORMAT 'PARQUET', CODEC 'ZSTD');
    "
    
    # Adding Transcript Metadata
    duckdb -c "    COPY (
                    SELECT 
                       geneDB.*,
                        tbl_transcript.Start AS Transcript_Start,
                        tbl_transcript.Stop AS Transcript_Stop,
                        tbl_transcript.Exon_count AS Exon_count
                        FROM read_parquet(${transcript_metadata}) AS tbl_transcript
                        RIGHT JOIN read_parquet('tmp_gene_constraints.parquet') AS geneDB
                        USING (Transcript_ID)
                ) TO "tmp_transcript_metadata.parquet" (FORMAT 'PARQUET', CODEC 'ZSTD');
    "

    
    # Adding Overlap Region
    duckdb -c " COPY (
                        SELECT geneDB.*, tbl_region_overlap.* EXCLUDE(CNV_ID),
                        -- Extract Chr from CNV_ID by taking substring before first '_'
                        SUBSTR(CNV_ID, 1, INSTR(CNV_ID || '_', '_') - 1) AS Chr
                        FROM read_csv_auto('${region_overlap}', delim='\\t', header=true) AS tbl_region_overlap
                        RIGHT JOIN read_parquet('tmp_transcript_metadata.parquet') AS geneDB
                        USING (CNV_ID)
                ) TO "geneDB.parquet" (FORMAT 'PARQUET', CODEC 'ZSTD', KV_METADATA {genome_version : \'${genome_version}\'});
    "
    """
}




workflow VEP_ANNOTATE {
    take:
    cnv_ch
    genome_version
    genomic_regions
    vep_cache
    gnomad_sv
    gnomad_constraints


    main:

    identifyUniqCNV(cnv_ch)

    if(genome_version == "GRCh38"){
        tbl_metadata = Channel.fromPath("${projectDir}/resources/Transcript_Metadata/transcriptDB_GRCh38.parquet")
        VEP_GRCh38(identifyUniqCNV.out, vep_cache, file(gnomad_sv) )
        vep_ch = VEP_GRCh38.out.results
    } else if(genome_version == "GRCh37") {
        tbl_metadata = Channel.fromPath("${projectDir}/resources/Transcript_Metadata/transcriptDB_GRCh37.parquet")
        VEP_GRCh37(identifyUniqCNV.out, vep_cache, file(gnomad_sv))
        vep_ch = VEP_GRCh37.out.results
    }

    computeOverlapRegion(identifyUniqCNV.out, genome_version, genomic_regions)

    db = buildGeneDB(vep_ch, gnomad_constraints, tbl_metadata, computeOverlapRegion.out, genome_version)
    emit:
    db
}