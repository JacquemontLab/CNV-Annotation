# Recurrent CNV Annotation (geneset_per_rCNV.tsv)

**Date:** 28/05/2026

This dataset of Regions of Interest (ROIs) for recurrent CNVs (rCNVs) was created primarily using [Clinical Genome Resource](https://www.clinicalgenome.org/), including all documented recurrent CNVs, with the addition of specific genes of interest. Some regions were manually curated; see the column **Reference** for details.

For each rCNV, we identify a geneset depending on the genome version (GRCh37 or GRCh38). To do so, we annotated the rCNV regions with **Ensembl VEP (Variant Effect Predictor)** to identify which transcripts fall within its boundary.


### Gene filtering Criteria

After VEP annotation, we kept only transcripts and genes meeting the following conditions for geneset construction:

* **Canonical transcripts** only
* **Protein-coding genes**
* **Variants/transcripts overlapping >50%** with the CNV region (OverlapPC > 50)
* **Transcripts not overlapping problematic regions** as defined by UCSC (ProblematicRegions_Overlap <50%)


### Recurrent CNVs identification

Two methods are used (see section 'Recurrent CNVs identification' of the main README for details):

* **Region method** (`bin/forge_rCNV.py`, run first): for the large recurrent regions only (`Locus_Group = 'Reccurent_CNV'`), an rCNV is flagged when the union of a sample's CNVs of the same Type covers >= 70% of the region coordinates (`GRCh37` / `GRCh38` columns). If no single CNV reaches 70% alone, a forged CNV (min Start, max End) is added.
* **Geneset method** (`bin/annotate_rCNV.py`, has priority): a CNV is flagged has recurrent if it overlaps all the genes in the geneset of a given rCNV_ID (considering only canonical transcripts of protein-coding genes) from resources/rCNV/geneset_per_rCNV.tsv .
For a given rCNV, its geneset is constructed based on the protein-coding canonical transcripts that it overlaps at 50% (see resources/rCNV/README.md for details). If more than one rCNV_ID is identified for a given CNV, then only the one with the largest geneset is kept (ties: alphabetical `rCNV_ID`).

Note: the `Locus_Group` value is spelled `Reccurent_CNV` in the file and the region method filters on this exact spelling. Some coordinates contain stray spaces (e.g. `X: 48447780-52444264`, `4:337779 -2009235`); they are removed when parsed.

`report_rCNV_method.pdf` is a short report validating this method based on the frequency of rCNVs observed in the literature and in a legacy approach previously used in the lab.


#### Scripts used

- rCNV_in_bed.R : Transforms rCNV coordinates from docs/recurrent_CNV_dataset.xlsx to a BED format.
- identify_gene_of_interest.sh : This script annotate the rCNV coordinates, keeping only genes and transcript of interest.
- identify_geneset.R : From the previous annotation construct geneset per rCNV, and save it.


#### Supplementary

- /docs/frequency_comparison.xlsx contains rCNV frequency (in percentage) compiled by Cecile Poulain, from different analysis/papers.
* Cécile Poulain’s response to reviewer (labelled: **supp\_cecile\_2025**)
* Kendall et al., 2019: [PMC6520248, Supplementary Table 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC6520248/#sec4) (labelled: **Kendall et al., 2019**)
* Crawford et al., 2019: [JMG, Supplementary Table 4](https://jmg.bmj.com/content/56/3/131.long) (labelled: **Crawford et al., 2019**)
* Legacy script in the lab (labelled: **historic\_lab**)
