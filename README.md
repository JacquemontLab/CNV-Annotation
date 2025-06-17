# CNV-DB-Builder

Nextflow pipeline for building a database from a single CNV file. The input that is expected is a bed-like format with at least the following columns:

```mermaid 
erDiagram
    direction TB
    CNV-Input{
        string SampleID
        string Chr
        int Start
        int End 
        string TYPE
    }
```

Where TYPE is a string that is either "DEL" or "DUP". Header names are optional while positioning is not. All other columns are passed over to the output.

### Dependencies 
 - python 3.13+
 - polars 
 - duckdb 
 - vep 113
 - Nextflow 25.04.2 


### Output
Minimally the output table is as follows:

| __dTYPE__ | __Column__ | __Description__                                    | 
|:--------- | -----------| -------------------------------------------------- |
|string     | CNV_ID             | ID of the CNV in the format of 'CHR_Start_End_TYPE'|
|string     | SampleID           | Cohort Specific ID for individual samples          |
|string     | Chr                | Chromosome Id. Optionally prefixed with 'Chr'      |
|int        | Start              | Chromosome start position. Ideally coordinates should match ensembl in that they are one-based and inclusive.|
|int        | End                | Chromosome End position.
|string     | TYPE               | CNV type. Either __'DEL'__ or __'DUP'__                    | 
|...| *__INPUT COLUMNS__* |                           |
|float      | segmentaldup_Overlap | Percentage base-pair overlap between CNV and segmental duplication regions. |
|string     | Gene_ID             | Ensembl ID for the overlapping gene with the CNV. |
|string     | Transcript_ID       | Ensembl ID for the __transcript__ overlapping with the CNV |
|string[]   | Consequence         | String list of Gene disruptions annotated by VEP.   | 
|boolean    | CANONICAL           | Transcript level canonical flag.                 |
|string     | MANE                | Matched Annotation from NCBI and EMBL-EBI (MANE) flag. [https://www.ncbi.nlm.nih.gov/refseq/MANE/](https://www.ncbi.nlm.nih.gov/refseq/MANE/). __Only available in HG38__ |
|string     | EXON                | String representation of the exons impacted by the CNV in the format of "<start_exon>-<end_exon>/<exon_count>" | 
|string     | INTRON              | String representation of the introns impacted by the CNV formatted as "<start_intron>-<end_intron>/<intron_count>" |
|float      | Exon_Overlap        | Percentage transform of the EXON column. See notes |
|float      | Transcript_bp_Overlap | Base-pair percentage overlap of the transcript with the CNV. |
|float      | Gnomad_Max_AF         | Maximum allele frequency of matching structural variant across populations. See notes. |  
|float 	    | LOEUF		    | From gnomAD V4:upper bound of 90% confidence interval for o/e ratio for high confidence pLoF variants (lower values indicate more constrained)|	

All other columns from the input are passed with their types estimated by python polars. 
The pipeline currently produces two parquet files, one being the formatted input file and the other the formatted VEP output,  and merges them based on CNV_ID alone. The structure between the ID columns are multiplcative and hierchical. SampleIDs, for instance will be duplicated in the following manner: 

```
# of dup SampleIDs = 1 * #CNV * #GENE * #Transcript
```   


### Notes
#### MANE 
MANE flag for transcript. Only supported in Hg38.
#### Gnomad_Max_AF 

Gnomad Allele Frequency (AF) annotations  for structural variants (SVs) are specific to the genome version.

__Hg19__ uses Gnomad V2 SV sites from here:
 https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz
    
- The fields extracted from the file are as follows:
    - AFR_AF
    - AMR_AF
    - EAS_AF
    - EUR_AF 

 __Hg38__ uses Gnomad V4 SV sites derived from WGS. The file was downloaded from here: https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz
 
 - The fields extracted from the file are as follows:
    - AF_nfe
    - AF_afr
    - AF_amr
    - AF_fin
    - AF_sas
    - AF_eas
    - AF_asj


A 70% reciprocal alignment is required for the CNV to be matched with a known SV. The maximum frequency is taken across all populations. In the event multiple gnomad SV annotations match, the maximum allele frequency is taken across SVs.

#### Exon_Overlap

By default, VEP reports CNVs that overlap with an exon in this format

    "<first_exon> - <last_exon> / <total_exon_count>"



Where "2-3/4" is a CNV that overlaps from the second to the third exon in gene of 4 exons. In order to convert this to a percentage format we apply the following function:

    Exon_Overlap = (<last_exon> - <first_exon> + 1) / <total_exon_count>

#### Transcript_bp_Overlap

This is a default field supplied by VEP. It is simply the base pair overlap the CNV shares with a transcript.




