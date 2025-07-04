#!/usr/bin/env nextflow 



//Formatting inputs for VEP annotation.
//This expects output from DigCNV, TO BE CHANGED 
process formatVEPInput {
    label 'quick'
    label 'polars'
    
    input:
    path cnvs 

    output:
    path "uniq_cnvs.bed"

    script:
    """
    prepare_cnvs_vep.py ${cnvs} "uniq_cnvs.bed"
    """
}



//Annotating for Gene disruptions
//process parameters in config. Requires apptainer image
process VEP_38 {

    label 'vep'
    input:
    path uniq_cnvs
    path vep_cache
    path gnomad_sv


    output:
    path "vep_out.tsv", emit : results
    path "*html", emit : summary
    path "*_warnings.txt", emit : warnings
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


process VEP_37 {
    label 'vep'
    
    input:
    path uniq_cnvs
    path vep_cache
    path gnomad_sv


    output:
    path "vep_out.tsv", emit : results
    path "*html", emit : summary
    path "*_warnings.txt", emit : warnings
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

process buildGeneDB {
    label 'quick'
    label 'polars'

    input:
    path vep_out
    path constraints
    path transcript_coors
    val genome_version

    output:
    path "geneDB.parquet"

    script:
    """
    #First transform large VEP output to parquet
    duckdb -c "COPY (SELECT * FROM read_csv(${vep_out}, delim = '\\t')) 
               TO 'tmp_db.parquet' (FORMAT 'PARQUET', CODEC 'ZSTD');"

    #Formatting output
    gene_db.py tmp_db.parquet tmp_formatted.parquet

    #Adding constraints file via right join on geneDB using gene_IDs
    duckdb -c "COPY ( 
                        SELECT geneDB.*, CAST(NULLIF(\\"lof.oe_ci.upper\\", 'NA') AS DOUBLE) AS LOEUF
                        FROM read_csv(\\"${constraints}\\", delim = '\\t') AS gnomad
                        RIGHT JOIN (SELECT * FROM read_parquet('tmp_formatted.parquet')) AS geneDB ON geneDB.Transcript_ID = gnomad.transcript
                        
        ) TO "tmp_gene_constraints.parquet" (FORMAT 'PARQUET', CODEC 'ZSTD');
    "
    #Adding transcript Coordinates
    duckdb -c " COPY (
                    SELECT 
                       geneDB.*,
                        {'Chr': t_coors.Chr, 'Start': t_coors.Start, 'Stop': t_coors.Stop} AS Transcript_Coordinates
                        FROM read_parquet(${transcript_coors}) AS t_coors
                        RIGHT JOIN read_parquet('tmp_gene_constraints.parquet') AS geneDB
                        USING (Transcript_ID)
                ) TO "geneDB.parquet" (FORMAT 'PARQUET', CODEC 'ZSTD', KV_METADATA {genome_version : \'${genome_version}\'});
    
    "
    """

}




workflow VEP_ANNOTATE {
    take:
    cnv_ch
    genome_version
    vep_cache
    gnomad_sv
    gnomad_constraints


    main:

    formatVEPInput(cnv_ch)

    if(genome_version == "GRCh38"){
        t_coors = Channel.fromPath("${projectDir}/resources/transcript_coords/transcript_coords_38.parquet")
        VEP_38(formatVEPInput.out, vep_cache, file(gnomad_sv) )
        vep_ch = VEP_38.out.results
    } else if(genome_version == "GRCh37") {
        t_coors = Channel.fromPath("${projectDir}/resources/transcript_coords/transcript_coords_37.parquet")
        VEP_37(formatVEPInput.out, vep_cache, file(gnomad_sv))
        vep_ch = VEP_37.out.results
    }


    db = buildGeneDB(vep_ch, gnomad_constraints, t_coors, genome_version)
    emit:
    db
    //vep_comments = vep_ch.out.comments
    //vep_summary  = vep_ch.out.summary

    
}
