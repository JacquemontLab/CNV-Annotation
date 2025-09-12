#!/usr/bin/env nextflow

/*
CNV_DB_Builder Nextflow Workflow
--------------------------------
This workflow performs the following steps:

1. Identify unique CNVs to reduce redundant queries for VEP annotation.
2. Compute overlap of CNVs with genomic regions.
3. Build a CNV database (Parquet format) combining CNV data with region annotations.
4. Annotate CNVs using VEP (Variant Effect Predictor) and generate LOEUF reports.
5. Produce summary PDFs for CNV and gene data.
6. Generate a run summary including duration, input, and output info.

Requirements:
- Nextflow DSL2
- Python scripts: prepare_cnvs_vep.py, add_regions_overlap.sh, format_overlap.sh, merge_cnv_with_region.py, pdf_dictionnary.py
- Polars library for Python
- VEP cache directory
*/

nextflow.enable.dsl=2
nextflow.preview.output = true
nextflow.enable.moduleBinaries = true

// Get Git hash at workflow launch
params.git_hash = "git -C ${projectDir} rev-parse HEAD".execute().text.trim()

// Default VEP dir following install script
params.vep_cache = "${projectDir}/resources"
params.genomic_regions = "${projectDir}/resources/Genome_Regions/Genome_Regions_data.tsv"
params.recurrent_path = "${projectDir}/resources/rCNV/geneset_per_rCNV.tsv"
params.gnomad_dir = "${params.vep_cache}/homo_sapiens" 

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
include { LOEUF_REPORT } from './modules/loeuf_report'
include { RCNV_ANNOTATION } from './modules/rCNV_annotation'


// It extracts unique CNV coordinates to reduce redundant queries
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



// Compute overlap of CNVs with genomic regions and add CNV_ID
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


// Merge CNV dataset with region overlaps to build CNV database
process buildCnvDB {
    label 'quick'

    input:
    path cnvs
    path region_overlap

    output:
    path "cnvDB_region.parquet"

    script:
    """
    # Adding Overlap Region
    merge_cnv_with_region.py ${cnvs} ${region_overlap} cnvDB_region.parquet
    """
}


// Generate summary PDFs from Parquet files
process produceSummaryPDF {
    label 'polars_duckdb'

    memory {
        // Dynamically allocate memory based on system total
        def memKB = new File('/proc/meminfo')
            .readLines()
            .find { it.startsWith('MemTotal:') }
            .replaceAll(/\D+/, '') as long
        return (memKB / 1024 / 1024).toInteger().GB
    }
    input:
    path parquet_input

    output:
    path "*_dictionary.pdf"

    script:
    """
    pdf_dictionary.py ${parquet_input} ${task.cpus} ${task.memory}
    """
}


// Build a launch summary file with workflow metadata and timing
process buildSummary {
    label 'quick'
    
    input:
    val cohort_tag
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

        // Step 2: Compute overlaps of CNVs with genomic regions
        computeOverlapRegion(uniq_cnv_ch, params.genome_version, params.genomic_regions)

        // Step 3: Merge CNVs with overlap information into a CNV database (Parquet format)
        buildCnvDB(cnvs_ch, computeOverlapRegion.out)

        // Step 4: Annotate CNVs using VEP (Variant Effect Predictor)
        VEP_ANNOTATE(
            uniq_cnv_ch,
            params.genome_version,
            params.vep_cache, 
            gnomad_AF,
            gnomad_constraints
        )

        // Step 5: Generate LOEUF-related figure using CNV DB and VEP annotation results
        LOEUF_REPORT(
            gnomad_constraints, //loeuf_metadata
            buildCnvDB.out,     //cnvDB
            VEP_ANNOTATE.out    //geneDB
        )
        
        RCNV_ANNOTATION(
            buildCnvDB.out,
            VEP_ANNOTATE.out,
            params.recurrent_path,
            params.genome_version)

        // Step 6: Produce PDF reports for CNV and gene annotation results
        pdf_cnv_ch = producePDFWorkflowCNV(RCNV_ANNOTATION.out.cnvDB_rCNV)
        pdf_gene_ch = producePDFWorkflowGene(VEP_ANNOTATE.out)
        
        // Step 7: Build a general summary report for the workflow run
        buildSummary(
            params.cohort_tag,
            params.cnvs,
            params.genome_version,
            params.git_hash,
            pdf_cnv_ch
        )

    // --- Publish outputs ---
    publish:
        cnv_db       = RCNV_ANNOTATION.out.cnvDB_rCNV          // Final CNV database
        gene_db      = VEP_ANNOTATE.out        // Annotated gene database
        summary      = buildSummary.out        // General workflow summary
        pdf_cnv      = pdf_cnv_ch              // CNV PDF report
        pdf_gene     = pdf_gene_ch             // Gene annotation PDF report
        loeuf_figure = LOEUF_REPORT.out        // LOEUF figures
}


// --- Output Sections ---
output {
    cnv_db {
        mode 'copy'
        path "${params.cohort_tag}/"
    }

    gene_db {
        mode 'copy'
        path "${params.cohort_tag}/"
    }

    pdf_gene {
        mode 'copy'
        path "${params.cohort_tag}/docs"
    }

    pdf_cnv {
        mode 'copy'
        path "${params.cohort_tag}/docs"
    }

    summary {
        mode 'copy'
        path "${params.cohort_tag}/docs/"
    }

    loeuf_figure {
        mode 'copy'
        path "${params.cohort_tag}/docs/"
    }
}
