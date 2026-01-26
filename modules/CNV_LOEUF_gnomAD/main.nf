#!/usr/bin/env nextflow


// This process computes, for each CNV_ID, the sum of LOEUF values
// from canonical transcripts (CANONICAL = TRUE) overlapping exons (Exon_Overlap > 0),
// then merge this LOEUF summary with the CNV database and add gnomAD Max AF.
process sum_loeuf_gnomAD_cnv {
    
    input:
    path uniq_cnv_parquet 
    path geneDB 
    path cnv_gnomAD

    output:
    path 'cnvDB_Loeuf_with_gnomAD.parquet'

    script:
    """
    duckdb -c "
        COPY (
            SELECT
                c.*,
                g.sum_LOEUF
            FROM '${uniq_cnv_parquet}' AS c
            LEFT JOIN (
                SELECT
                    CNV_ID,
                    SUM(LOEUF) AS sum_LOEUF
                FROM '${geneDB}'
                WHERE CANONICAL = TRUE
                AND Exon_Overlap > 0
                GROUP BY CNV_ID
            ) AS g
            USING (CNV_ID)
        ) TO 'cnvDB_Loeuf.parquet' (FORMAT 'PARQUET');
        "

    duckdb -c "
        COPY (
            SELECT
                cnvDB.*,
                gnomad.Gnomad_Max_AF
            FROM read_parquet('cnvDB_Loeuf.parquet') AS cnvDB
            LEFT JOIN read_parquet('${cnv_gnomAD}') AS gnomad
            USING (CNV_ID)
        )
        TO 'cnvDB_Loeuf_with_gnomAD.parquet' (FORMAT PARQUET);
    "
    """
}