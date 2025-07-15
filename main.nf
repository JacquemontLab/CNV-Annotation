#!/usr/bin/env nextflow

nextflow.preview.output = true
nextflow.enable.moduleBinaries = true

def gnomad_AF
def gnomad_constraints = "${params.vep_cache}/ressources_LOEUF/gnomad.v4.1.constraint_metrics.tsv"


switch (params.genome_version) {
    case "GRCh38":
        gnomad_AF = "${params.vep_cache}/ressources_gnomAD/gnomad.v4.1.sv.sites.vcf.bgz"
        break
    case "GRCh37":
        gnomad_AF = "${params.vep_cache}/ressources_gnomAD/gnomad.v2.1.sv.sites.vcf.bgz"
        break
    default:
        error "Unsupported genome version '${params.genome_version}'. Use 'GRCh38' or 'GRCh37'."
}


include { VEP_ANNOTATE } from './modules/vep_annotate'


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
}

process joinTables {
    tag { "${tag}" }

    input:
    path vep_db
    path input_db
    val tag

    output:
    path "cnvDB.parquet"

    script:
    """
    # Extract chromosomes to process into chr_list.txt
    duckdb -c "COPY (SELECT DISTINCT Chr FROM '${input_db}') TO 'chr_list.txt' WITH (FORMAT csv, HEADER false);"

    # Read chromosomes into an array
    mapfile -t CHR_LIST < chr_list.txt

    echo "Chromosomes to process: \${CHR_LIST[*]}"

    for CHR in "\${CHR_LIST[@]}"; do
        echo "Processing chromosome \$CHR"
        
        # Extract available memory (GiB) as integer
        export AVAILABLE_MEM=\$(free -h | awk '/Mem:/ {gsub(/Gi/, "", \$7); print int(\$7)}')
        LIMIT_MEM=\$((AVAILABLE_MEM - 20))

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

        COPY cnvDB_all TO 'cnvDB.parquet' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 500000,
                                KV_METADATA {DB_Run_Name: \'${workflow.runName}\'});
    "

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
    val genome_version
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
       genome_version: ${genome_version}
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
    log.info "Using genome version: ${params.genome_version}"
    log.info "gnomAD AF file: ${gnomad_AF}"
    log.info "gnomAD constraint file: ${gnomad_constraints}"

    main:
        //gnomad    = file( projectDir / "resources" / params.genome_version / "gnomad.v*sv.sites.vcf.bgz")
    
        cnvs_ch   = Channel.fromPath(file(params.cnvs))

        tag_file_size = cnvs_ch.map { file ->
            def sizeGB = file.size() / 1_073_741_824
            sizeGB > 2 ? 'LargeFile' : 'SmallFile'
        }

        VEP_ANNOTATE (  cnvs_ch, 
                        params.genome_version,
                        params.genomic_regions, 
                        params.vep_cache, 
                        gnomad_AF,
                        gnomad_constraints )
        
        buildCnvDB    ( cnvs_ch)

        joinTables    ( VEP_ANNOTATE.out,
                        buildCnvDB.out,
                        tag_file_size)
        
        produceSummaryPDF (
                        joinTables.out,
                        "64",
                        "3.5" )
        
        buildSummary  (params.cohort_tag,
                        params.cnvs,
                        params.genome_version,
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