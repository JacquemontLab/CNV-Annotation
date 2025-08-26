
# Getting Transcript Metadata from GTF file:
```
curl https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz > Homo_sapiens.GRCh38.113.gtf.gz

duckdb -c "
    CREATE TABLE tbl AS (
        SELECT * FROM read_csv('Homo_sapiens.GRCh38.113.gtf.gz', delim = '\t', all_varchar = true)
    );

    -- Add Transcript_ID by extracting it from column8
    ALTER TABLE tbl ADD COLUMN Transcript_ID VARCHAR;
    UPDATE tbl 
    SET Transcript_ID = REGEXP_EXTRACT(column8, 'transcript_id \"(.*?)\"', 1);


    -- Add Gene_ID by extracting it from column8
    ALTER TABLE tbl ADD COLUMN Gene_ID VARCHAR;
    UPDATE tbl 
    SET Gene_ID = REGEXP_EXTRACT(column8, 'gene_id \"(.*?)\"', 1);

    -- Extract transcript_biotype
    ALTER TABLE tbl ADD COLUMN Transcript_biotype VARCHAR;
    UPDATE tbl 
    SET Transcript_biotype = REGEXP_EXTRACT(column8, 'transcript_biotype \"(.*?)\"', 1);

    -- Extract transcript_source
    ALTER TABLE tbl ADD COLUMN Transcript_source VARCHAR;
    UPDATE tbl 
    SET Transcript_source = REGEXP_EXTRACT(column8, 'transcript_source \"(.*?)\"', 1);

    -- Count exons per Transcript_ID
    CREATE TABLE exon_counts AS
    SELECT Transcript_ID, COUNT(*) AS exon_count
    FROM tbl
    WHERE column2 = 'exon'
    GROUP BY Transcript_ID;

    -- Export transcript entries with exon counts
    COPY (
        SELECT 
            column0::VARCHAR AS Chr,
            column3::INTEGER AS Start,
            column4::INTEGER AS Stop,
            tbl.Transcript_ID,
            tbl.Gene_ID,
            COALESCE(exon_counts.exon_count, 0) AS Exon_count,
            CASE WHEN Transcript_biotype = 'protein_coding' THEN 1 ELSE 0 END AS Is_coding,
            tbl.Transcript_source
        FROM tbl
        LEFT JOIN exon_counts USING (Transcript_ID)
        WHERE column2 = 'transcript'
    ) TO 'transcriptDB_GRCh38.parquet';
"


curl https://ftp.ensembl.org/pub/grch37/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz > Homo_sapiens.GRCh37.87.gtf.gz


duckdb -c "
    CREATE TABLE tbl AS (
        SELECT * FROM read_csv('Homo_sapiens.GRCh37.87.gtf.gz', delim = '\t', all_varchar = true)
    );

    -- Add Transcript_ID by extracting it from column8
    ALTER TABLE tbl ADD COLUMN Transcript_ID VARCHAR;
    UPDATE tbl 
    SET Transcript_ID = REGEXP_EXTRACT(column8, 'transcript_id \"(.*?)\"', 1);


    -- Add Gene_ID by extracting it from column8
    ALTER TABLE tbl ADD COLUMN Gene_ID VARCHAR;
    UPDATE tbl 
    SET Gene_ID = REGEXP_EXTRACT(column8, 'gene_id \"(.*?)\"', 1);

    -- Extract transcript_biotype
    ALTER TABLE tbl ADD COLUMN Transcript_biotype VARCHAR;
    UPDATE tbl 
    SET Transcript_biotype = REGEXP_EXTRACT(column8, 'transcript_biotype \"(.*?)\"', 1);

    -- Extract transcript_source
    ALTER TABLE tbl ADD COLUMN Transcript_source VARCHAR;
    UPDATE tbl 
    SET Transcript_source = REGEXP_EXTRACT(column8, 'transcript_source \"(.*?)\"', 1);

    -- Count exons per Transcript_ID
    CREATE TABLE exon_counts AS
    SELECT Transcript_ID, COUNT(*) AS exon_count
    FROM tbl
    WHERE column2 = 'exon'
    GROUP BY Transcript_ID;

    -- Export transcript entries with exon counts
    COPY (
        SELECT 
            column0::VARCHAR AS Chr,
            column3::INTEGER AS Start,
            column4::INTEGER AS Stop,
            tbl.Transcript_ID,
            tbl.Gene_ID,
            COALESCE(exon_counts.exon_count, 0) AS Exon_count,
            CASE WHEN Transcript_biotype = 'protein_coding' THEN 1 ELSE 0 END AS Is_coding,
            tbl.Transcript_source
        FROM tbl
        LEFT JOIN exon_counts USING (Transcript_ID)
        WHERE column2 = 'transcript'
    ) TO 'transcriptDB_GRCh37.parquet';
"


```



