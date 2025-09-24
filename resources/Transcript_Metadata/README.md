
# Getting Transcript Metadata from GTF file:

```bash
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

    -- Add gene_name by extracting it from column8
    ALTER TABLE tbl ADD COLUMN Gene_Name VARCHAR;
    UPDATE tbl 
    SET Gene_Name = REGEXP_EXTRACT(column8, 'gene_name \"(.*?)\"', 1);

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
            tbl.Gene_ID,
            tbl.Gene_Name,
            tbl.Transcript_ID,
            tbl.Transcript_biotype,
            tbl.Transcript_source,
            COALESCE(exon_counts.exon_count, 0) AS Exon_count
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

    -- Add gene_name by extracting it from column8
    ALTER TABLE tbl ADD COLUMN Gene_Name VARCHAR;
    UPDATE tbl 
    SET Gene_Name = REGEXP_EXTRACT(column8, 'gene_name \"(.*?)\"', 1);

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
            tbl.Gene_ID,
            tbl.Gene_Name,
            tbl.Transcript_ID,
            tbl.Transcript_biotype,
            tbl.Transcript_source,
            COALESCE(exon_counts.exon_count, 0) AS Exon_count
        FROM tbl
        LEFT JOIN exon_counts USING (Transcript_ID)
        WHERE column2 = 'transcript'
    ) TO 'transcriptDB_GRCh37.parquet';
"

```

### Getting Transcript Overlaps With Problematic Regions

#### GRCh38


```bash
cat ../Genome_Regions/Genome_Regions_data.tsv | grep problematic_regions | grep GRCh38 | cut -f1-3 > pRegions_38.bed
cat pRegions_38.bed | tr --delete "chr"   > pRegions_38_rename.bed
bedtools intersect -b pRegions_38_rename.bed -a Homo_sapiens.GRCh38.113.gtf -wao  > pRegion38.gtf.bed
```
Extract transcript entries and reduce to final output
```sql
CREATE OR REPLACE TABLE inter AS (SELECT * FROM read_csv('pRegion38.gtf.bed', delim = '\t', all_varchar = true) WHERE column02 == 'transcript');
ALTER TABLE inter ADD COLUMN Transcript_ID VARCHAR;
ALTER TABLE inter ADD COLUMN Transcript_problematic_regions_Overlap FLOAT;
UPDATE inter SET Transcript_ID = regexp_extract(column08,  'transcript_id \"(.*?)\"',  1);
UPDATE inter SET Transcript_problematic_regions_Overlap =  (CAST(column12 AS INTEGER) / (CAST(column04 AS INTEGER) - CAST(column03 AS INTEGER) + 1));

.shell cp transcriptDB_GRCh38.parquet transcriptDB_GRCh38_BACKUP.parquet

COPY (
    SELECT tDB.*, 
        inter_max.Transcript_problematic_regions_Overlap,
    FROM ( 
        SELECT
            Transcript_ID,
            MAX(Transcript_problematic_regions_Overlap) AS Transcript_problematic_regions_Overlap
                FROM inter 
                GROUP BY Transcript_ID 
            ) AS inter_max
    RIGHT JOIN 'transcriptDB_GRCh38_BACKUP.parquet' AS tDB 
    USING (Transcript_ID)) TO 'transcriptDB_GRCh38.parquet';



```

#### GRCh37
```bash
curl https://ftp.ensembl.org/pub/grch37/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz > Homo_sapiens.GRCh37.87.gtf.gz
cat ../Genome_Regions/Genome_Regions_data.tsv | grep problematic_regions | grep GRCh37  | cut -f1-3 > pRegions_37.bed
cat pRegions_37.bed | tr --delete "chr"   > pRegions_37_rename.bed
bedtools intersect -b pRegions_37_rename.bed -a Homo_sapiens.GRCh37.87.gtf.gz -wao  > pRegion37.gtf.bed
```

```sql
CREATE OR REPLACE TABLE inter AS (SELECT * FROM read_csv('pRegion37.gtf.bed', delim = '\t', all_varchar = true) WHERE column02 == 'transcript');
ALTER TABLE inter ADD COLUMN Transcript_ID VARCHAR;
ALTER TABLE inter ADD COLUMN Transcript_problematic_regions_Overlap FLOAT;
UPDATE inter SET Transcript_ID = regexp_extract(column08,  'transcript_id \"(.*?)\"',  1);
UPDATE inter SET Transcript_problematic_regions_Overlap =  (CAST(column12 AS INTEGER) / (CAST(column04 AS INTEGER) - CAST(column03 AS INTEGER) + 1));

.shell cp transcriptDB_GRCh37.parquet transcriptDB_GRCh37_BACKUP.parquet
COPY (
    SELECT tDB.*, 
        inter_max.Transcript_problematic_regions_Overlap,
    FROM ( 
        SELECT
            Transcript_ID,
            MAX(Transcript_problematic_regions_Overlap) AS Transcript_problematic_regions_Overlap 
                FROM inter 
                GROUP BY Transcript_ID 
            ) AS inter_max
    RIGHT JOIN 'transcriptDB_GRCh37_BACKUP.parquet' AS tDB 
    USING (Transcript_ID)) TO 'transcriptDB_GRCh37.parquet';

```