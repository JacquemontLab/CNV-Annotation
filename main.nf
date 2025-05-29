#!/usr/bin/env nextflow

nextflow.preview.output = true

//test params
params.cnvs           = "${projectDir}/test-data/CNVs.tsv"
params.genome_version = "GRCh37"
params.cohort_tag     = "ALSPAC"

//TODO replace CNVdb builder process with a polars script for adding cnv_id and TYPE flag 

include { VEP_ANNOTATE } from './modules/vep_annotate'


/* process buildSampleDB {
    label 'quick'

    input:
    path plink
    path qc_scores


    output:
    path "sampleDB.parquet"

    script:
    """
    duckdb -c "COPY (
                    SELECT *
                    FROM  read_csv(${plink}, delim='\\t', header=True) AS t1
                    RIGHT JOIN read_csv(${qc_scores}, delim='\\t', header=True) AS t2
                    USING (SampleID)
                    ) 
                TO 'sampleDB.parquet' (FORMAT 'parquet');"
    """

    stub:
    """
    touch sampleDB.parquet
    """
} */

process buildCnvDB {
    
    label 'quick'

    input:
    path cnvs
    val genome_version
    path regions_file 

    output:
    path "cnvDB.parquet"

    script:
    """
    add_regions_overlap.sh ${cnvs} ${regions_file} ${genome_version} "CNVs_with_genomic_regions_overlap.tsv"

    cnv_db_builder_lite.py CNVs_with_genomic_regions_overlap.tsv cnvDB.parquet
    """

    stub:
    """
    touch cnvDB.parquet
    """
}

process buildSummary {
    label 'quick'
    input:
    val cohort_tag

    output:
    path "launch_report.txt"

    script:
    """
       cat <<EOF > launch_report.txt
       CNV_DB_Builder ${cohort_tag} run summary:
       run name: ${workflow.runName}
       version: ${workflow.manifest.version}
       configs: ${workflow.configFiles}
       workDir: ${workflow.workDir}
       launch_user: ${workflow.userName}
       start_time: ${workflow.start}
       duration: ${workflow.duration}

       Command:
       ${workflow.commandLine}

    
    """

    stub:
    """
    touch launch_report.txt
    """
}
workflow {

    main:
        //gnomad    = file( projectDir / "resources" / params.genome_version / "gnomad.v*sv.sites.vcf.bgz")
    
        cohort_ch = params.cohort_tag
        cnvs_ch   = Channel.fromPath(file(params.cnvs))



        VEP_ANNOTATE (  cnvs_ch, 
                        params.genome_version, 
                        params.vep_cache, 
                        params.gnomad_AF )


       /*  buildSampleDB ( plink_ch,
                        qc_ch ) */
        
        buildCnvDB    ( cnvs_ch,
                        params.genome_version,
                        params.genome_regions,
                                                )
        buildSummary  (params.cohort_tag)

    publish:
        cnv_gene_db  = VEP_ANNOTATE.out
        //sample_db    = buildSampleDB.out
        cnv_db       = buildCnvDB.out
        summary      = buildSummary.out


}

output {
    cnv_gene_db {
        mode 'copy'
        path "${params.cohort_tag}/"
    }
/*     sample_db {
        path "${params.cohort_tag}/"
    } */
    cnv_db {
        mode 'copy'
        path "${params.cohort_tag}/"
    }
    summary {
        mode 'copy'
        path "${params.cohort_tag}/"
    }
    
}

