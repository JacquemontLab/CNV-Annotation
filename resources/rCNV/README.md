# Recurrent CNV Annotation

**Date:** 26/08/2025

This dataset of Regions of Interest (ROIs) for recurrent CNVs (rCNVs) was created primarily using [Clinical Genome Resource](https://www.clinicalgenome.org/), including all documented recurrent CNVs, with the addition of specific genes of interest. Some regions were manually curated; see the column **Reference** for details.

For each rCNV, we identify a geneset depending on the genome version (GRCh37 or GRCh38). To do so, we annotated the rCNV regions with **Ensembl VEP (Variant Effect Predictor)** to identify which transcripts fall within its boundary.


### Gene filtering Criteria

After VEP annotation, we kept only transcripts and genes meeting the following conditions for geneset construction:

* **Canonical transcripts** only
* **Protein-coding genes**
* **Variants/transcripts overlapping >50%** with the CNV region (OverlapPC > 50)
* **Transcripts not overlapping problematic regions** as defined by UCSC (overlap <50%)


#### Scripts used

- rCNV_in_bed.R : Transforms rCNV coordinates from docs/recurrent_CNV_dataset.xlsx to a BED format.
- identify_gene_of_interest.sh : This script annotate the rCNV coordinates, keeping only genes and transcript of interest.
- identify_geneset.R : From the previous annotation construct geneset per rCNV, and save it.


#### Supplementary

- /docs/frequency_comparison.xlsx contains rCNV frequency compiled by Cecile Poulain, from different analysis/papers.