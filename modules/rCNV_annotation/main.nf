#!/usr/bin/env nextflow

// --- Process: annotate_rCNV ---
// This process annotates CNVs with gene information and recurrent CNV flags.
// Inputs:
//   - cnvDB: path to the CNV database (Parquet format)
//   - geneDB: path to the gene annotation database (Parquet format)
//   - recurrent_path: path to a TSV file containing recurrent CNV gene sets
//   - genome_version: genome build to use (e.g., GRCh37 or GRCh38)
// Outputs:
//   - cnvDB_rCNV.parquet: CNV database annotated with flagged recurrent CNVs
process annotate_rCNV {
    label 'polars_duckdb'

    input:
    path cnvDB
    path geneDB
    path recurrent_path
    val genome_version

    output:
    path 'cnvDB_rCNV.parquet', emit : cnvDB_rCNV

    script:
    """
    annotate_rCNV.py \
        --geneDB_path ${geneDB} \
        --cnvDB_path ${cnvDB} \
        --recurrent_path ${recurrent_path} \
        --cnvDB_flagged_parquet cnvDB_rCNV.parquet \
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

    emit:
    cnvDB_rCNV = results.cnvDB_rCNV
}
