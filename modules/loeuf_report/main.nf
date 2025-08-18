#!/usr/bin/env nextflow


// This process merges the CNV database with the Gene database using CNV_ID as the key.
// The output is a merged Parquet file containing both CNV and gene information.
process merge_cnv_gene {
    label 'quick'
    
    input:
    path cnvDB 
    path geneDB 

    output:
    path 'mergedDB.parquet'

    script:
    """
    duckdb -c "
        COPY (
            SELECT cnv.*,
                gene.*
            FROM read_parquet('${cnvDB}') AS cnv
            LEFT JOIN read_parquet('${geneDB}') AS gene
            USING (CNV_ID)
        ) TO 'mergedDB.parquet' (FORMAT 'PARQUET', CODEC 'ZSTD');
        "

    """
}

// Generates a LOEUF-based figure (CNV enrichment per LOEUF decile) from the merged CNV-Gene database.
process loeuf_report {
    label 'quick'
    
    input:
    path loeuf_metadata 
    path mergeDB 

    output:
    path "loeuf_report.png", emit : figure

    script:
    """
    loeuf_cnv_duckdb.py -c ${mergeDB} -l ${loeuf_metadata} -o loeuf_report.png -f Two_Algorithm_Overlap
    """
}


// --- Workflow: LOEUF_REPORT ---
// Main workflow to compute LOEUF report.
// Steps:
// 1. Merge CNV and Gene databases.
// 2. Generate LOEUF figure from the merged database.
workflow LOEUF_REPORT {
    take:
    loeuf_metadata
    cnvDB
    geneDB

    main:

    loeuf_report(loeuf_metadata, merge_cnv_gene(cnvDB, geneDB))

    emit:
    loeuf_report_png = loeuf_report.out.figure
}
