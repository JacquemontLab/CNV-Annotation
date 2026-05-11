#!/usr/bin/env nextflow


/*
CNV-Annotation Nextflow Pipeline
================================

Description
-----------
Builds an annotated CNV (copy-number variant) database from input CNVs. 
Reduces redundancy by identifying unique CNVs, annotates them with VEP, 
computes overlaps with genomic regions, aggregates LOEUF scores, integrates 
recurrent CNVs, and produces CNV- and gene-level PDF reports.

Workflow Steps
--------------
1. Identify unique CNVs to minimize redundant annotation.
2. Annotate CNVs using VEP (functional and gnomAD metrics).
3. Compute overlaps with genomic regions.
4. Aggregate LOEUF scores per CNV.
5. Build the CNV database (Parquet format).
6. Add recurrent CNV annotations.
7. Generate PDF reports and a global workflow summary.

Usage
-----
nextflow run main.nf --cnvs path/to/cnvs.tsv --regions path/to/regions.bed \
  --genome_version GRCh38 --vep_cache /path/to/vep_cache --outdir results
*/

nextflow.enable.dsl = 2

// Get Git hash at workflow launch
params.git_hash = "git -C ${projectDir} rev-parse HEAD".execute().text.trim()

// Default VEP dir following install script
params.vep_cache = "${projectDir}/resources/vep_cache"
params.genomic_regions = "${projectDir}/resources/Genome_Regions/Genome_Regions_data.tsv"
params.recurrent_path = "${projectDir}/resources/rCNV/geneset_per_rCNV.tsv"
params.recurrent_freq_path = "${projectDir}/resources/rCNV/docs/frequency_comparison.xlsx"

def gnomad_AF
def gnomad_constraints = "${params.vep_cache}/ressources_LOEUF/gnomad.v4.1.constraint_metrics.tsv"


// Select appropriate gnomAD file based on genome version
switch (params.genome_version) {
    case "GRCh38":
        gnomad_AF = "${params.vep_cache}/ressources_gnomAD/gnomad.v4.1.sv.sites.vcf.bgz"
        break
    case "GRCh37":
        gnomad_AF = "${params.vep_cache}/ressources_gnomAD/gnomad_v2.1_sv.sites.vcf.bgz" 
        break
    default:
        error "Unsupported genome version '${params.genome_version}'. Use 'GRCh38' or 'GRCh37'."
}


// Include external modules for VEP annotation and LOEUF report generation
include { VEP_ANNOTATE } from './modules/vep_annotate'
include { sum_loeuf_gnomAD_cnv } from './modules/CNV_LOEUF_gnomAD'
include { RCNV_ANNOTATION } from './modules/rCNV_annotation'


// It extracts unique CNV coordinates to reduce redundant queries
process identifyUniqCNV {
    
    input:
    path cnvs 

    output:
    path "uniq_cnvs.bed", emit: uniq_cnv_bed
    path "uniq_cnvs.parquet", emit: uniq_cnv_parquet

    script:
    """
    prepare_cnvs_vep.py ${cnvs} "uniq_cnvs.bed" "uniq_cnvs.parquet"
    """
}



// Compute overlap of CNVs with genomic regions and add CNV_ID
process computeOverlapRegion {
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


// Merge CNV dataset with region overlaps to build CNV database
process buildCnvDB {

    input:
    path cnvs
    path region_overlap
    path sum_loeuf_gnomAD_cnv
    path cnvDB_rCNV

    output:
    path "cnvDB.parquet"

    script:
    """
    # Step 1: Build CNV regions with overlap
    merge_cnv_with_region.py ${cnvs} ${region_overlap} cnvDB_region.parquet

    # Step 2: Merge with sum_loeuf_gnomAD_cnv database on CNV_ID
    duckdb -c "
        COPY (
        SELECT r.*, c.sum_LOEUF, c.Gnomad_Max_AF
        FROM read_parquet('cnvDB_region.parquet') AS r
        LEFT JOIN read_parquet('${sum_loeuf_gnomAD_cnv}') AS c
        USING (CNV_ID)
        )
        TO 'cnvDB_tmp.parquet' (FORMAT PARQUET);
    "

    # Step 3: Merge with rCNV database on CNV_ID
    duckdb -c "
        COPY (
        SELECT r.*, c.rCNV_ID
        FROM read_parquet('cnvDB_tmp.parquet') AS r
        LEFT JOIN read_parquet('${cnvDB_rCNV}') AS c
        USING (CNV_ID)
        )
        TO 'cnvDB.parquet' (FORMAT PARQUET);
    "
    """
}


// Generate summary PDFs from Parquet files
process produceSummaryPDF {

    input:
    path parquet_input

    output:
    path "*_columns_report.pdf"

    script:
    // Compute memory in GB if task.memory exists; else leave null
    def memGb = task.memory ? (task.memory.toGiga() * 0.85).intValue() : null

    """
    # If memGb > 0 (HPC), use it; else detect 85% of VM RAM
    if [ ${memGb} -gt 0 ]; then
        MEM_GB=${memGb}
    else
        MEM_GB=\$(free -k | awk '/^Mem:/ {print int(\$2*0.85/1024/1024)}')
    fi

    echo "Using memory limit: \${MEM_GB} GB"

    pdf_columns_report.py ${parquet_input} ${task.cpus} \${MEM_GB}
    """
}


process QCReportPDF {
    label "Rmarkdown"

    input:
    path cnv_dataset_report_rmd
    val dataset_name
    path cnvDB
    path geneDB
    path recurrent_path
    path recurrent_freq_path

    output:
    path "cnv_dataset_qc.pdf"

    script:
    """
    # The markdown is initially designed to run by providing paths to resource scripts from Git repositories
    # That's what we are reconstructing here below
    mkdir -p resources/rCNV/
    cp ${recurrent_path} resources/rCNV/

    mkdir -p resources/rCNV/docs/
    cp ${recurrent_freq_path} resources/rCNV/docs/
  
    mkdir -p bin/make_report/
    cp ${projectDir}/bin/cnv_trio_inheritance.py bin/make_report/

    Rscript -e "rmarkdown::render('${cnv_dataset_report_rmd}', 
        params=list(path_dataset='.',
        dataset_name='${dataset_name}',
        path_CNVANNOTATION_repo='.',
        path_CNVCALLER_repo='.'), 
        output_file='cnv_dataset_qc.pdf')"
    """
}


// Build a launch summary file with workflow metadata and timing
process buildSummary {
    input:
    val dataset_name
    val cnvs_path
    val genome_version
    val git_hash
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

    # Convert duration to hours, minutes, seconds
    hours=\$(( duration / 3600 ))
    minutes=\$(( (duration % 3600) / 60 ))
    seconds=\$(( duration % 60 ))

    cat <<EOF > launch_report.txt
    CNV-Annotation ${dataset_name} run summary:
    run name: ${workflow.runName}
    version: ${workflow.manifest.version}
    configs: ${workflow.configFiles}
    workDir: ${workflow.workDir}
    input_file: ${cnvs_path}
    genome_version: ${genome_version}
    launch_user: ${workflow.userName}
    start_time: ${workflow.start}
    duration: \${hours}h \${minutes}m \${seconds}s

    Command:
    ${workflow.commandLine}

    Git hash working version:
    commit ${git_hash}
    """

    stub:
    """
    touch launch_report.txt
    """
}



// --- Sub-workflows for producing PDFs ---
workflow producePDFWorkflowCNV {
    take:
        input_ch

    main:
        pdf_ch = produceSummaryPDF(input_ch)

    emit:
        pdf_ch
}

workflow producePDFWorkflowGene {
    take:
        input_ch

    main:
        pdf_ch = produceSummaryPDF(input_ch)

    emit:
        pdf_ch
}



// --- Main Workflow ---
// This is the core pipeline execution section, orchestrating the CNV DB building, annotation, and reporting.

workflow {
    // Log basic configuration
    log.info "Using genome version: ${params.genome_version}"
    log.info "gnomAD AF file: ${gnomad_AF}"
    log.info "gnomAD constraint file: ${gnomad_constraints}"

    main:    
        // Load the input CNV file(s) into a Nextflow channel
        cnvs_ch = Channel.fromPath(file(params.cnvs))

        // Step 1: Identify unique CNVs to reduce redundancy before annotation
        uniq_cnv_ch = identifyUniqCNV(cnvs_ch)

        // Step 2: Annotate CNVs using VEP (Variant Effect Predictor)
        VEP_ANNOTATE(
            uniq_cnv_ch.uniq_cnv_bed,
            params.genome_version,
            params.vep_cache, 
            gnomad_AF,
            gnomad_constraints
        )

        // Step 3: Compute overlaps of CNVs with genomic regions
        computeOverlapRegion(uniq_cnv_ch.uniq_cnv_bed, params.genome_version, params.genomic_regions)

        // Step 4: Compute sum of LOEUF values for CNVs and add gnomAD Max AF       
        sum_loeuf_gnomAD_cnv(
            uniq_cnv_ch.uniq_cnv_parquet,
            VEP_ANNOTATE.out.geneDB,
            VEP_ANNOTATE.out.cnv_gnomAD)

        // Step 5: Annotate CNVs with recurrent CNV information
        RCNV_ANNOTATION(
            uniq_cnv_ch.uniq_cnv_parquet,
            VEP_ANNOTATE.out.geneDB,
            params.recurrent_path,
            params.genome_version)

        // Step 6: Merge CNVs with overlap information into a CNV database (Parquet format)
        buildCnvDB(cnvs_ch, computeOverlapRegion.out, sum_loeuf_gnomAD_cnv.out, RCNV_ANNOTATION.out.cnvDB_rCNV)

        // Step 7: Produce PDF reports for CNV and gene annotation results
        pdf_cnvDB = producePDFWorkflowCNV(buildCnvDB.out)
        pdf_geneDB = producePDFWorkflowGene(VEP_ANNOTATE.out.geneDB)
        
        // Step 8: Produce PDF  QC Report of the CNVs Dataset
        QCReportPDF(
            file("${projectDir}/bin/cnv_dataset_report.Rmd"),
            params.dataset_name,
            buildCnvDB.out,
            VEP_ANNOTATE.out.geneDB,
            params.recurrent_path,
            params.recurrent_freq_path)

        // Step 9: Build a general summary report for the workflow run
        buildSummary(
            params.dataset_name,
            params.cnvs,
            params.genome_version,
            params.git_hash,
            QCReportPDF.out
        )

    // --- Publish outputs ---
    publish:
        cnv_db        = buildCnvDB.out                     // Final CNV database
        gene_db       = VEP_ANNOTATE.out.geneDB            // Annotated gene database
        summary       = buildSummary.out                   // General workflow summary
        pdf_cnv       = pdf_cnvDB                          // CNV PDF report
        pdf_gene      = pdf_geneDB                         // Gene annotation PDF report
        pdf_qc_report = QCReportPDF.out                    // QC Report – CNVs Dataset
}


// --- Output Sections ---
output {
    cnv_db {
        mode 'copy'
        path "${params.dataset_name}/"
    }

    gene_db {
        mode 'copy'
        path "${params.dataset_name}/"
    }

    pdf_gene {
        mode 'copy'
        path "${params.dataset_name}/docs"
    }

    pdf_cnv {
        mode 'copy'
        path "${params.dataset_name}/docs"
    }

    pdf_qc_report {
        mode 'copy'
        path "${params.dataset_name}/docs"
    }

    summary {
        mode 'copy'
        path "${params.dataset_name}/docs/"
    }
}
