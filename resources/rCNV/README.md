# Recurrent CNV Annotation

**Date:** 26/08/2025

This dataset of Regions of Interest (ROIs) for recurrent CNVs (rCNVs) was created primarily using [Clinical Genome Resource](https://www.clinicalgenome.org/), including all documented recurrent CNVs, with the addition of specific genes of interest. Some regions were manually curated; see the column **Reference** for details.

For each rCNV, we identify a geneset depending on the genome version (GRCh37 or GRCh38). To do so, we annotated the rCNV regions with **Ensembl VEP (Variant Effect Predictor)** to identify which transcripts fall within each CNV boundary.

### Filtering Criteria

After annotation, we retained only transcripts and genes meeting the following conditions:

* **Canonical transcripts** only
* **Protein-coding genes**
* **Variants/transcripts overlapping >50%** with the CNV region (OverlapPC > 50)
* **Transcripts not overlapping problematic regions** as defined by UCSC (overlap <50%)
