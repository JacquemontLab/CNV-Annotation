#!/usr/bin/env python3
"""
forge_rCNV.py

Purpose:
    Region-based recurrent CNV (rCNV) annotation, run at the very start of the
    pipeline on the raw cohort CNV file. Fragmented calls covering an rCNV region
    are merged into a single "forged" CNV that then goes through the whole pipeline.

Method (per SampleID x Type x rCNV region):
    1. Overlap every CNV with the rCNV regions (coordinates of the chosen genome version),
       restricted to Locus_Group = 'Reccurent_CNV' (large recurrent regions, not single genes).
    2. Coverage = union of the sample's CNVs (same Type) clipped to the region / region length.
    3. A region passes if coverage >= --min_coverage (default 0.7).
    4. If several passing regions overlap each other for the same SampleID x Type,
       only the largest region is kept.
    5. For each kept region:
         - If ONE CNV alone covers >= min_coverage of the region: no forging, that CNV
           is annotated (rCNV_ID, rCNV_ID_method = 'region', rCNV_forged = FALSE).
         - Otherwise a forged CNV is added: Start = min(Start), End = max(End) of all
           CNVs overlapping the region; every other column is copied from the largest
           of those CNVs (rCNV_ID, rCNV_ID_method = 'region', rCNV_forged = TRUE).
           The original CNVs are kept untouched (rCNV_forged = FALSE, no rCNV_ID).

Output:
    Input TSV (all original columns, values unchanged) + forged rows, with 3 new columns:
    rCNV_ID (e.g. 22q11.21_del), rCNV_ID_method ('region' or empty), rCNV_forged (true/false).

Usage:
    forge_rCNV.py --cnvs cnvs.tsv --recurrent_path geneset_per_rCNV.tsv \
                  --genome_version GRCh38 --output cnvs_forged.tsv
"""

import argparse
import duckdb


def main(args):
    con = duckdb.connect(database=":memory:")

    # Case-insensitive lookup of the required columns (same convention as prepare_cnvs_vep.py)
    with open(args.cnvs) as f:
        header = f.readline().rstrip("\n").split("\t")
    col = {name.lower(): name for name in header}
    SAMPLE, CHR, START, END, TYPE = (f'"{col[k]}"' for k in ["sampleid", "chr", "start", "end", "type"])

    # 1. Load CNVs as text (values are written back unchanged) + a row id
    con.execute(f"""
    CREATE TABLE cnv AS
    SELECT row_number() OVER () AS _row_id, *
    FROM read_csv('{args.cnvs}', delim='\t', header=true, all_varchar=true);
    """)

    con.execute(f"""
    CREATE TABLE cnv_pos AS
    SELECT _row_id,
           {SAMPLE} AS sample,
           {TYPE}   AS type,
           regexp_replace({CHR}, '^chr', '') AS chrom,
           CAST({START} AS BIGINT) AS s,
           CAST({END}   AS BIGINT) AS e
    FROM cnv;
    """)

    # 2. Load rCNV regions, coordinates like '22:17392953-18591860' (spaces removed)
    #    Only the large recurrent CNV regions (Locus_Group = 'Reccurent_CNV', spelling as in the file);
    #    single-gene regions are left to the geneset method.
    con.execute(f"""
    CREATE TABLE regions AS
    SELECT rCNV_ID,
           regexp_replace(split_part(coord, ':', 1), '^chr', '') AS chrom,
           CAST(split_part(split_part(coord, ':', 2), '-', 1) AS BIGINT) AS r_start,
           CAST(split_part(split_part(coord, ':', 2), '-', 2) AS BIGINT) AS r_end,
           r_end - r_start AS region_len
    FROM (
        SELECT rCNV_ID, replace({args.genome_version}, ' ', '') AS coord
        FROM read_csv('{args.recurrent_path}', delim='\t', header=true, all_varchar=true)
        WHERE Locus_Group = 'Reccurent_CNV'
    );
    """)

    # 3. CNV x region overlaps, with the CNV clipped to the region (os, oe)
    con.execute("""
    CREATE TABLE ov AS
    SELECT c._row_id, c.sample, c.type, c.chrom, c.s, c.e,
           r.rCNV_ID, r.r_start, r.r_end, r.region_len,
           greatest(c.s, r.r_start) AS os,
           least(c.e, r.r_end)      AS oe
    FROM cnv_pos c
    JOIN regions r
      ON c.chrom = r.chrom AND c.s < r.r_end AND c.e > r.r_start;
    """)

    # 4. Union coverage per sample x type x region (merge overlapping clipped CNVs)
    con.execute("""
    CREATE TABLE coverage AS
    WITH ordered AS (
        SELECT *,
               max(oe) OVER (PARTITION BY sample, type, rCNV_ID ORDER BY os, oe
                             ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING) AS prev_max_oe
        FROM ov
    ),
    islands AS (
        SELECT *,
               sum(CASE WHEN prev_max_oe IS NULL OR os > prev_max_oe THEN 1 ELSE 0 END)
                   OVER (PARTITION BY sample, type, rCNV_ID ORDER BY os, oe
                         ROWS UNBOUNDED PRECEDING) AS island
        FROM ordered
    ),
    island_len AS (
        SELECT sample, type, rCNV_ID, island, max(oe) - min(os) AS len
        FROM islands
        GROUP BY sample, type, rCNV_ID, island
    )
    SELECT sample, type, rCNV_ID, sum(len) AS covered_bp
    FROM island_len
    GROUP BY sample, type, rCNV_ID;
    """)

    # 5. Passing candidates (union coverage >= threshold)
    con.execute(f"""
    CREATE TABLE candidates AS
    SELECT o.sample, o.type, o.rCNV_ID, o.chrom, o.r_start, o.r_end, o.region_len,
           cv.covered_bp / o.region_len           AS union_cov,
           max(o.oe - o.os) / o.region_len        AS best_single_cov,
           arg_max(o._row_id, o.oe - o.os)        AS best_single_row,  -- CNV covering most of the region
           arg_max(o._row_id, o.e - o.s)          AS largest_row,      -- largest CNV (source of other columns)
           min(o.s) AS forged_start,
           max(o.e) AS forged_end,
           count(*) AS n_cnv
    FROM ov o
    JOIN coverage cv USING (sample, type, rCNV_ID)
    GROUP BY o.sample, o.type, o.rCNV_ID, o.chrom, o.r_start, o.r_end, o.region_len, cv.covered_bp
    HAVING cv.covered_bp / o.region_len >= {args.min_coverage};
    """)

    # 6. Among overlapping passing regions (same sample x type), keep the largest one
    #    (ties broken alphabetically on rCNV_ID)
    con.execute(f"""
    CREATE TABLE selected AS
    SELECT c.*,
           c.rCNV_ID || '_' || lower(c.type) AS rCNV_ID_with_type,
           c.best_single_cov < {args.min_coverage} AS forged
    FROM candidates c
    WHERE NOT EXISTS (
        SELECT 1 FROM candidates o
        WHERE o.sample = c.sample AND o.type = c.type AND o.chrom = c.chrom
          AND o.rCNV_ID <> c.rCNV_ID
          AND o.r_start < c.r_end AND o.r_end > c.r_start
          AND (o.region_len > c.region_len
               OR (o.region_len = c.region_len AND o.rCNV_ID < c.rCNV_ID))
    );
    """)

    # 7a. Not forged: annotate the existing CNV (one region per CNV, the largest)
    con.execute("""
    CREATE TABLE annotated AS
    SELECT best_single_row AS _row_id, rCNV_ID_with_type AS rCNV_ID
    FROM selected
    WHERE NOT forged
    QUALIFY row_number() OVER (PARTITION BY best_single_row
                               ORDER BY region_len DESC, rCNV_ID) = 1;
    """)

    # 7b. Forged: copy the largest CNV, replace Start/End with min/max of the overlapping CNVs
    con.execute(f"""
    CREATE TABLE forged_rows AS
    SELECT c.* EXCLUDE (_row_id)
             REPLACE (CAST(f.forged_start AS VARCHAR) AS {START},
                      CAST(f.forged_end   AS VARCHAR) AS {END}),
           f.rCNV_ID_with_type AS rCNV_ID,
           'region'            AS rCNV_ID_method,
           true                AS rCNV_forged
    FROM selected f
    JOIN cnv c ON c._row_id = f.largest_row
    WHERE f.forged;
    """)

    # 8. Original CNVs (annotated or not, input order kept) + forged CNVs at the end
    con.execute(f"""
    COPY (
        SELECT * EXCLUDE (_order)
        FROM (
            SELECT c.* EXCLUDE (_row_id),
                   a.rCNV_ID,
                   CASE WHEN a.rCNV_ID IS NOT NULL THEN 'region' END AS rCNV_ID_method,
                   false AS rCNV_forged,
                   c._row_id AS _order
            FROM cnv c
            LEFT JOIN annotated a USING (_row_id)
            UNION ALL BY NAME
            SELECT *, NULL AS _order FROM forged_rows
        )
        ORDER BY _order NULLS LAST
    ) TO '{args.output}' (HEADER, DELIMITER '\t');
    """)

    # Log a short summary
    n_annot, n_forged = con.execute(
        "SELECT (SELECT count(*) FROM annotated), (SELECT count(*) FROM forged_rows)"
    ).fetchone()
    print(f"rCNV region method: {n_annot} CNVs annotated without forging, {n_forged} forged CNVs added.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Region-based rCNV annotation and CNV forging (DuckDB)")
    parser.add_argument("--cnvs", required=True, help="Input cohort CNV TSV (SampleID, Chr, Start, End, Type, ...)")
    parser.add_argument("--recurrent_path", required=True, help="rCNV regions TSV (geneset_per_rCNV.tsv)")
    parser.add_argument("--genome_version", required=True, choices=["GRCh37", "GRCh38"])
    parser.add_argument("--output", required=True, help="Output TSV with forged CNVs and rCNV columns")
    parser.add_argument("--min_coverage", type=float, default=0.7, help="Min fraction of the region covered (default 0.7)")
    main(parser.parse_args())
