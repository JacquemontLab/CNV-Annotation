#!/usr/bin/env nextflow

nextflow.preview.output = true
nextflow.enable.moduleBinaries = true

//test params
params.cnvs           = "${projectDir}/test-data/CNVs.tsv"
params.genome_version = "GRCh37"
params.cohort_tag     = "ALSPAC"



include { VEP_ANNOTATE } from './modules/vep_annotate'


process buildSampleDB {
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
}

process buildCnvDB {
    
    label 'quick'
    label 'polars'

    input:
    path cnvs
    val genome_version
    path regions_file 

    output:
    path "inputDB.parquet"

    script:
    """
    add_regions_overlap.sh ${cnvs} ${regions_file} ${genome_version} "CNVs_with_genomic_regions_overlap.tsv"

    cnv_db_builder_lite.py CNVs_with_genomic_regions_overlap.tsv inputDB.parquet
    """

    stub:
    """
    touch cnvDB.parquet
    """
}

process joinTables {
    label 'quick' // expect to run on launch job with a good number of cpus ~16 minimum and 32gb ram

    input:
    path vep_db
    path input_db

    output:
    path "cnvDB.parquet"

    script:
    """
    duckdb -c "
        COPY(
            SELECT * FROM ${input_db}
            FULL JOIN (SELECT * EXCLUDE (Allele, Location) FROM ${vep_db})
            USING (CNV_ID)
        ) TO "cnvDB.parquet" (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 5_000_000,
                              KV_METADATA {DB_Run_Name: \'${workflow.runName}\'});
    "
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
                        params.gnomad_AF,
                        params.gnomad_constraints )


       /*  buildSampleDB ( plink_ch,
                        qc_ch ) */
        
        buildCnvDB    ( cnvs_ch,
                        params.genome_version,
                        params.genome_regions,
                                                )

        joinTables    ( VEP_ANNOTATE.out,
                        buildCnvDB.out          )

        buildSummary  (params.cohort_tag)

    publish:
        cnv_db       = joinTables.out
        summary      = buildSummary.out


}

output {
    cnv_db {
        mode 'copy'
        path "${params.cohort_tag}/"
    }
    summary {
        mode 'copy'
        path "${params.cohort_tag}/"
    }
    
}

