
# Genomic regions of interest TSV File only for canonical chromosomes
# All the following data are combined in one unique file Genome_Regions_data.tsv
Chr\tStart\tEnd\tRegion\tGenomeVersion

## PAR (Pseudoautosomal Region) regions dataset
From https://www.ncbi.nlm.nih.gov/grc/human on 25/04/2025

GRCh37.p13
chrX	60001	2699520	PAR1	GRCh37
chrX	154931044	155260560	PAR2	GRCh37

GRCh38.p14
chrX	10001	2781479	PAR1	GRCh38
chrX	155701383	156030895	PAR2	GRCh38


## X-Transpose region dataset
From Timothy H Webster, Madeline Couse, Bruno M Grande, Eric Karlins, Tanya N Phung, Phillip A Richmond, Whitney Whitford, Melissa A Wilson, Identifying, understanding, and correcting technical artifacts on the sex chromosomes in next-generation sequencing data, GigaScience, Volume 8, Issue 7, July 2019, giz074, https://doi.org/10.1093/gigascience/giz074

"We define the XTR on the X chromosome as beginning at the start of DXS1217 and ending at the end of DXS3"
 
GRCh37.p13
https://grch37.ensembl.org/Homo_sapiens/Marker/Details?m=sWXD902
https://grch37.ensembl.org/Homo_sapiens/Marker/Details?m=sWXD298

DXS1217	chromosome X:88395845-88396079	GRCh37
DXS3	chromosome X:92582890-92583067	GRCh37

GRCh38.p14
https://useast.ensembl.org/Homo_sapiens/Marker/Details?db=core;m=DXS1217;r=X:89140845-89141079
https://useast.ensembl.org/Homo_sapiens/Marker/Details?db=core;m=DXS3;r=X:93327891-93328068

DXS1217	chromosome X:89140845-89141079	GRCh38
DXS3	chromosome X:93327891-93328068	GRCh38

### XTR region coordinates
chrX	88395845	92583067	XTR	GRCh37
chrX	89140845	93328068	XTR	GRCh38

## Problematic Region

Problematic_GRCh38.tsv from :
https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=2763451842_9nXejNOmv3oAIqSDNs99CqacdGPH&clade=mammal&org=Human&db=hg38&hgta_group=allTracks&hgta_track=problematic&hgta_table=comments&hgta_regionType=genome&position=chr7%3A155%2C799%2C529-155%2C812%2C871&hgta_outputType=primaryTable&hgta_outFileName=Problematic_GRCh38.tsv

Considering Tables UCSC Unusual Regions, ENCODE Blacklist V2, GRC Exclusions.

For GRCh37.tsv Tables : UCSC Unusual Regions, ENCODE Blacklist V2, all GIAB and all NCBI zone proposed.
https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=2763451842_9nXejNOmv3oAIqSDNs99CqacdGPH&clade=mammal&org=Human&db=hg19&hgta_group=allTracks&hgta_track=problematic&hgta_table=filterSSE&hgta_regionType=genome&position=chr7%3A155%2C592%2C223-155%2C605%2C565&hgta_outputType=primaryTable&hgta_outFileName=problematic_GRCh37.tsv


After downloading every tables :
```bash
for f in *38.tsv; do tail -n +2 "$f"; done | cut -f1-3 | awk -F'\t' -v OFS='\t' '{print $0, "problematic_regions", "GRCh38"}' > Problematic_GRCh38_regions.tsv

for f in *37.tsv; do tail -n +2 "$f"; done | cut -f1-3 | awk -F'\t' -v OFS='\t' '{print $0, "problematic_regions", "GRCh37"}' > Problematic_GRCh37_regions.tsv
```
