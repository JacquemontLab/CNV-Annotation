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

    input:
    path cnvs

    output:
    path "inputDB.parquet"

    script:
    """
    cnv_db_builder_lite.py ${cnvs} inputDB.parquet
    """

    stub:
    """
    touch cnvDB.parquet
    """
}

process joinTables {
    label 'high_memory' // expect to run on launch job with a good number of cpus ~16 minimum and 32gb ram

    input:
    path vep_db
    path input_db

    output:
    path "cnvDB.parquet"

    script:
    """
    # Extract chromosomes to process into chr_list.txt
    duckdb -c "COPY (SELECT DISTINCT Chr FROM 'inputDB.parquet') TO 'chr_list.txt' WITH (FORMAT csv, HEADER false);"

    # Read chromosomes into an array
    mapfile -t CHR_LIST < chr_list.txt

    echo "Chromosomes to process: \${CHR_LIST[*]}"

    for CHR in "\${CHR_LIST[@]}"; do
        echo "Processing chromosome \$CHR"
        
        # Extract available memory (GiB) as integer
        export AVAILABLE_MEM=\$(free -h | awk '/Mem:/ {gsub(/Gi/, "", \$7); print int(\$7)}')
        LIMIT_MEM=\$((AVAILABLE_MEM - 4))

        echo "Available memory: \$AVAILABLE_MEM GiB"
        echo "Memory limit after subtraction: \$LIMIT_MEM GiB"

        # Run DuckDB with memory limit
        duckdb -c " 
            SET memory_limit='\${LIMIT_MEM}GB';

            COPY(
                SELECT * 
                FROM (SELECT * FROM '${input_db}' WHERE Chr = '\${CHR}')
                FULL JOIN (
                    SELECT * EXCLUDE (Allele, Location, Chr) FROM '${vep_db}' WHERE Chr = '\${CHR}'
                )
                USING (CNV_ID)
            )  TO 'cnvDB_\${CHR}.parquet' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 500000,
                                KV_METADATA {DB_Run_Name: \'${workflow.runName}\'});
        "
        if [ \$? -ne 0 ]; then
            echo "DuckDB failed on chromosome \$CHR — likely due to OOM"
            exit 100  
        fi
    done

    # Merge parquet files into one (concatenate)
    duckdb cnvDB.duckdb -c "
        CREATE TABLE cnvDB_all AS SELECT * FROM read_parquet('cnvDB_chr1.parquet') LIMIT 0;
    "

    for file in cnvDB_*.parquet; do
        echo "Appending \$file ..."
        duckdb cnvDB.duckdb -c "
        INSERT INTO cnvDB_all SELECT * FROM read_parquet('\$file');
        "
    done

    duckdb cnvDB.duckdb -c "
        SET memory_limit='\${LIMIT_MEM}GB';

        COPY cnvDB_all TO 'cnvDB.parquet' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 500000);
    "

    if [ \$? -ne 0 ]; then
        echo "DuckDB failed during merging parquet files"
        exit 101
    fi

    """

    stub:
    """
    touch cnvDB.parquet
    """
}


process produceSummaryPDF {

    input:
    path cnvDB_parquet
    val cpu
    val mem_per_cpu

    output:
    path "cnvDB_dictionary.pdf"

    script:
    """
    pdf_dictionnary.py ${cnvDB_parquet} ${cpu} ${mem_per_cpu}
    """
}


process buildSummary {
    label 'quick'
    input:
    val cohort_tag
    val cnvs_path
    path last_outfile

    output:
    path "launch_report.txt"

    script:
    """
        # Convert workflow start datetime to epoch seconds
        start_sec=\$(date -d "${workflow.start}" +%s)
        # Get current time in epoch seconds
        end_sec=\$(date +%s)

        # Calculate duration in seconds
        duration=\$(( end_sec - start_sec ))

        # Convert duration to minutes and seconds
        minutes=\$(( duration / 60 ))
        seconds=\$(( duration % 60 ))

       cat <<EOF > launch_report.txt
       CNV_DB_Builder ${cohort_tag} run summary:
       run name: ${workflow.runName}
       version: ${workflow.manifest.version}
       configs: ${workflow.configFiles}
       workDir: ${workflow.workDir}
       input_file: ${cnvs_path}
       launch_user: ${workflow.userName}
       start_time: ${workflow.start}
       duration: \${minutes} minutes and \${seconds} seconds

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
    
        cnvs_ch   = Channel.fromPath(file(params.cnvs))

        VEP_ANNOTATE (  cnvs_ch, 
                        params.genome_version,
                        params.genomic_regions, 
                        params.vep_cache, 
                        params.gnomad_AF,
                        params.gnomad_constraints )
        
        buildCnvDB    ( cnvs_ch)

        joinTables    ( VEP_ANNOTATE.out,
                        buildCnvDB.out)
        
        produceSummaryPDF (
                        joinTables.out,
                        "64",
                        "3.5" )
        
        buildSummary  (params.cohort_tag,
                        params.cnvs,
                        produceSummaryPDF.out )

    publish:
        cnv_db       = joinTables.out
        summary      = buildSummary.out
        pdf_summary  = produceSummaryPDF.out


}

output {
    cnv_db {
        mode 'copy'
        path "${params.cohort_tag}/"
    }

    pdf_summary {
        mode 'copy'
        path "${params.cohort_tag}/docs/"
    }

    summary {
        mode 'copy'
        path "${params.cohort_tag}/docs/"
    }
    
}