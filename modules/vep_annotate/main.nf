#!/usr/bin/env nextflow

// ================================================================
// CNV Annotation with VEP (Variant Effect Predictor)
// ---------------------------------------------------------------
// This workflow annotates unique CNVs against genes using VEP.
// It supports both GRCh37 and GRCh38 genome assemblies.
// It integrates gnomAD SV data, transcript metadata, and LOEUF
// constraints to produce a gene-level Parquet database.
// ================================================================


// ---------------------------
// Process: VEP_GRCh38
// ---------------------------
// Runs VEP for GRCh38 assembly, generates tabular and HTML outputs,
// and extracts comment lines for logs.
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

    # detect CPUs inside the container
    CPUS=\$(nproc)
    echo "Using \$CPUS CPUs for VEP"

    vep -i ${uniq_cnvs} -o vep_out.tsv\
    -cache\
    --tab\
    --dir_cache ${vep_cache}\
    --offline\
    --force_overwrite\
    --numbers\
    --fork \$CPUS \
    --biotype\
    --overlaps\
    --canonical\
    --mane\
    --max_sv_size 100000000\
    --verbose\
    --assembly GRCh38 \
    --custom file="./${gnomad_sv}",short_name=gnomad,format=VCF,reciprocal=1,overlap_cutoff=70,same_type=1,fields=AF_nfe%AF_afr%AF_amr%AF_fin%AF_sas%AF_eas%AF_asj \
    --fields "Uploaded_variation,Location,Allele,Gene,Feature,Consequence,BIOTYPE,CANONICAL,MANE,EXON,INTRON,OverlapPC,gnomad_AF_nfe,gnomad_AF_afr,gnomad_AF_amr,gnomad_AF_fin,gnomad_AF_sas,gnomad_AF_eas,gnomad_AF_asj"


    grep -E '^\\s*#' vep_out.tsv > vep_comments.txt
    """
}


// ---------------------------
// Process: VEP_GRCh37
// ---------------------------
// Runs VEP for GRCh37 assembly with genome-specific gnomAD fields.
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
    
    # detect CPUs inside the container
    CPUS=\$(nproc)
    echo "Using \$CPUS CPUs for VEP"

    vep -i ${uniq_cnvs} -o vep_out.tsv\
    -cache\
    --tab\
    --dir_cache ${vep_cache}\
    --offline\
    --force_overwrite\
    --numbers\
    --fork \$CPUS \
    --biotype\
    --overlaps\
    --canonical\
    --max_sv_size 100000000\
    --verbose\
    --assembly GRCh37 \
    --custom file="./${gnomad_sv}",short_name=gnomad,format=VCF,reciprocal=1,overlap_cutoff=70,same_type=1,fields=AFR_AF%AMR_AF%EAS_AF%EUR_AF \
    --fields "Uploaded_variation,Location,Allele,Gene,Feature,Consequence,BIOTYPE,CANONICAL,MANE,EXON,INTRON,OverlapPC,gnomad_AFR_AF,gnomad_AMR_AF,gnomad_EAS_AF,gnomad_EUR_AF"


    grep -E '^\\s*#' vep_out.tsv > vep_comments.txt
    """
}


// Integrates VEP output, gnomAD constraints (LOEUF), and transcript metadata.
// Produces a compressed Parquet file containing genome-version metadata.
process buildGeneDB {

    input:
    path vep_out
    path gnomad_constraints
    path transcript_metadata
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
                        tbl_transcript.Gene_Name AS Gene_Name,
                        tbl_transcript.Start AS Transcript_Start,
                        tbl_transcript.Stop AS Transcript_Stop,
                        tbl_transcript.Exon_count AS Exon_count
                        FROM read_parquet(${transcript_metadata}) AS tbl_transcript
                        RIGHT JOIN read_parquet('tmp_gene_constraints.parquet') AS geneDB
                        USING (Transcript_ID)
                ) TO "geneDB.parquet" (FORMAT 'PARQUET', CODEC 'ZSTD');
    "

    """
}




// ---------------------------
// Workflow: VEP_ANNOTATE
// ---------------------------
// Chooses genome assembly-specific VEP process, then builds gene database.
workflow VEP_ANNOTATE {
    take:
    uniq_cnvs
    genome_version
    vep_cache
    gnomad_sv
    gnomad_constraints


    main:

    if(genome_version == "GRCh38"){
        transcript_metadata = Channel.fromPath("${projectDir}/resources/Transcript_Metadata/transcriptDB_GRCh38.parquet")
        VEP_GRCh38(uniq_cnvs, vep_cache, file(gnomad_sv) )
        vep_ch = VEP_GRCh38.out.results
        
    } else if(genome_version == "GRCh37") {
        transcript_metadata = Channel.fromPath("${projectDir}/resources/Transcript_Metadata/transcriptDB_GRCh37.parquet")
        VEP_GRCh37(uniq_cnvs, vep_cache, file(gnomad_sv))
        vep_ch = VEP_GRCh37.out.results
    }

    db = buildGeneDB(vep_ch, gnomad_constraints, transcript_metadata, genome_version)

    emit:
    db
}
