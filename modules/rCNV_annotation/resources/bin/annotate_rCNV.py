#!/usr/bin/env python3
 
import duckdb
import argparse
import os

# Connect to DuckDB (in-memory)
con = duckdb.connect(database=':memory:')
    
    
def create_table_from_file(table_name, file_path):
    ext = os.path.splitext(file_path)[1].lower()
    if ext in ['.tsv', '.csv']:
        # read CSV/TSV
        con.execute(f"""
            CREATE TABLE {table_name} AS
            SELECT *
            FROM read_csv_auto('{file_path}');
        """)
    elif ext in ['.parquet', '.parq']:
        # read Parquet
        con.execute(f"""
            CREATE TABLE {table_name} AS
            SELECT *
            FROM read_parquet('{file_path}');
        """)
    else:
        raise ValueError(f"Unsupported file extension: {ext}")



def main(args):

    # 1. Load input
    create_table_from_file("geneDB", args.geneDB_path)
    create_table_from_file("cnvDB", args.cnvDB_path)
    create_table_from_file("recurrent", args.recurrent_path)

    # 2. Explode geneDB
    con.execute("""
    CREATE TABLE geneDB_exploded AS
    SELECT CNV_ID, Gene_ID, Allele AS Type
    FROM geneDB
    WHERE Exon_Overlap > 0
      AND CANONICAL = 'true';
    """)

    # 3. Explode recurrent
    gene_col = f"geneset_{args.genome_version}"  # e.g., geneset_GRCh38 or geneset_GRCh37

    con.execute(f"""
    CREATE TABLE recurrent_exploded AS
    SELECT 
        rCNV_ID,
        TRIM(g.Gene_ID) AS Gene_ID
    FROM recurrent
    CROSS JOIN UNNEST(string_split({gene_col}, ',')) AS g(Gene_ID);
    """)

    # 4. Filter
    con.execute("""
    CREATE TABLE geneDB_filtered AS
    SELECT g.*
    FROM geneDB_exploded g
    INNER JOIN recurrent_exploded r
        ON g.Gene_ID = r.Gene_ID;
    """)

    # 5. Matching counts
    con.execute("""
    CREATE TABLE matching_counts AS
    WITH cnv_gene_counts AS (
        SELECT 
            CNV_ID, 
            MIN(Type) AS Type,
            COUNT(DISTINCT Gene_ID) AS total_genes_per_CNV
        FROM geneDB_filtered
        GROUP BY CNV_ID
    )
    SELECT 
        g.CNV_ID,
        r.rCNV_ID,
        COUNT(DISTINCT g.Gene_ID) AS matched_genes,
        c.Type,
        c.total_genes_per_CNV
    FROM geneDB_filtered g
    JOIN recurrent_exploded r
        ON g.Gene_ID = r.Gene_ID
    JOIN cnv_gene_counts c
        ON g.CNV_ID = c.CNV_ID
    GROUP BY g.CNV_ID, r.rCNV_ID, c.Type, c.total_genes_per_CNV;
    """)

    # 6. Recurrent counts
    con.execute("""
    CREATE TABLE recurrent_counts AS
    SELECT rCNV_ID, COUNT(DISTINCT Gene_ID) AS total_genes
    FROM recurrent_exploded
    GROUP BY rCNV_ID;
    """)

    # 7. Full matches
    con.execute("""
    CREATE TABLE full_matches AS
    WITH ranked_matches AS (
        SELECT 
            m.CNV_ID,
            m.Type,
            m.rCNV_ID,
            m.matched_genes,
            m.total_genes_per_CNV,
            r.total_genes AS recurrent_total_genes,
            ROW_NUMBER() OVER (
                PARTITION BY m.CNV_ID
                ORDER BY m.matched_genes DESC
            ) AS rn
        FROM matching_counts m
        JOIN recurrent_counts r
          ON m.rCNV_ID = r.rCNV_ID
        WHERE m.matched_genes = r.total_genes
    )
    SELECT *
    FROM ranked_matches
    WHERE rn = 1;
    """)

    # 8. Add rCNV_ID with type suffix
    con.execute("""
    CREATE TABLE full_matches_with_rCNV AS
    SELECT 
        f.CNV_ID,
        f.rCNV_ID,
        r.rCNV_ID || '_' || LOWER(f.Type) AS rCNV_ID_with_type,
        f.matched_genes,
        f.total_genes_per_CNV,
        f.recurrent_total_genes
    FROM full_matches f
    JOIN recurrent r
        ON f.rCNV_ID = r.rCNV_ID;
    """)

    # 9. Join to cnvDB
    con.execute("""
    CREATE TABLE cnvDB_flagged AS
    SELECT g.*, r.rCNV_ID_with_type AS rCNV_ID
    FROM cnvDB g
    LEFT JOIN full_matches_with_rCNV r
      ON g.CNV_ID = r.CNV_ID;
    """)

    # 10. Save flagged cnvDB
    con.execute(f"""
    COPY cnvDB_flagged
    TO '{args.cnvDB_flagged_parquet}'
    (FORMAT PARQUET);
    """)


    # Expand recurrent into recurrent_expanded (with _dup and _del)
    con.execute("""
    CREATE OR REPLACE TABLE recurrent_expanded AS
    SELECT rCNV_ID || '_dup' AS rCNV_ID
    FROM recurrent
    UNION ALL
    SELECT rCNV_ID || '_del' AS rCNV_ID
    FROM recurrent;
    """)

    # Now join and count
    con.execute("""
    CREATE OR REPLACE TABLE sample_rCNV_counts AS
    SELECT 
        r.rCNV_ID,
        COALESCE(COUNT(DISTINCT c.SampleID), 0) AS num_samples
    FROM recurrent_expanded r
    LEFT JOIN cnvDB_flagged c
        ON c.rCNV_ID = r.rCNV_ID
    GROUP BY r.rCNV_ID
    ORDER BY r.rCNV_ID;
    """)
    
    # Save as TSV file
    con.execute(f"""
    COPY sample_rCNV_counts 
    TO '{args.recurrent_sample_counts}' 
    (HEADER, DELIMITER '\t');
    """)
    
    print("Processing complete!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DuckDB CNV-Recurrent processing")
    parser.add_argument("--geneDB_path", required=True, help="Input geneDB file (TSV or Parquet)")
    parser.add_argument("--cnvDB_path", required=True, help="Input cnvDB file (TSV or Parquet)")
    parser.add_argument("--recurrent_path", required=True, help="Input recurrent CNV gene set file (TSV)")
    parser.add_argument("--cnvDB_flagged_parquet", required=True, help="Output path for flagged cnvDB Parquet file")
    parser.add_argument("--recurrent_sample_counts", required=True, help="Output path for recurrent sample counts TSV")
    parser.add_argument("--genome_version", required=True, choices=["GRCh37", "GRCh38"], help="Genome version to use")
    args = parser.parse_args()

    
    main(args)
