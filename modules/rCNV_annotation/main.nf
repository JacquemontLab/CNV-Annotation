#!/usr/bin/env nextflow

// --- Process: annotate_rCNV ---
// This process annotates CNVs with gene information and recurrent CNV flags.
// Inputs:
//   - cnvDB: path to the CNV database (Parquet format)
//   - geneDB: path to the gene annotation database (Parquet format)
//   - recurrent_path: path to a TSV file containing recurrent CNV gene sets
//   - genome_version: genome build to use (e.g., GRCh37 or GRCh38)
// Outputs:
//   - cnvDB.parquet: CNV database annotated with flagged recurrent CNVs
//   - rCNV_sample_counts.tsv: table of sample counts per recurrent CNV
process annotate_rCNV {
    label 'quick'

    input:
    path cnvDB
    path geneDB
    path recurrent_path
    val genome_version

    output:
    path 'cnvDB.parquet', emit : cnvDB_rCNV
    path 'rCNV_sample_counts.tsv', emit : rCNV_sample_counts

    script:
    """
    annotate_rCNV.py \
        --geneDB_path ${geneDB} \
        --cnvDB_path ${cnvDB} \
        --recurrent_path ${recurrent_path} \
        --cnvDB_flagged_parquet cnvDB.parquet \
        --recurrent_sample_counts rCNV_sample_counts.tsv \
        --genome_version ${genome_version}
    """
}


// --- Workflow: RCNV_ANNOTATION ---
// Main workflow to annotate CNVs using gene and recurrent CNV information.
// Steps:
// 1. Call the annotate_rCNV process with the given CNV DB, gene DB, and recurrent file.
// 2. Emit a flaggedDB map containing both outputs for downstream use.
workflow RCNV_ANNOTATION {
    take:
    cnvDB
    geneDB
    recurrent_path
    genome_version

    main:
    // Call the process; returns a map of emitted outputs
    results = annotate_rCNV(cnvDB, geneDB, recurrent_path, genome_version)

    // Assign each emitted output to a variable
    cnvDB_rCNV = results.cnvDB_rCNV
    rCNV_sample_counts = results.rCNV_sample_counts

    emit:
    cnvDB_rCNV
    rCNV_sample_counts
}
