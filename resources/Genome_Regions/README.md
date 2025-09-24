
# Genomic regions of interest TSV File only for canonical chromosomes
# All the following data are combined in one unique file Genome_Regions_data.tsv
Chr\tStart\tEnd\tRegion\tGenomeVersion

## Segmental Duplication regions dataset

Downloaded on 25/04/2025 from https://genome.ucsc.edu/cgi-bin/hgTables

SegmentalDups_GRCh37.bed from :
https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=2529510726_A566RApY6cLEYg7x3NCXAq93ABrZ&clade=mammal&org=Human&db=hg19&hgta_group=allTracks&hgta_track=genomicSuperDups&hgta_table=0&hgta_regionType=genome&position=chr7%3A155%2C592%2C223-155%2C605%2C565&hgta_outputType=bed&hgta_outFileName=SegmentalDups_GRCh37.bed

SegmentalDups_GRCh38.bed from :
https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=2529510726_A566RApY6cLEYg7x3NCXAq93ABrZ&clade=mammal&org=&db=hg38&hgta_group=allTracks&hgta_track=genomicSuperDups&hgta_table=genomicSuperDups&hgta_regionType=genome&position=&hgta_outputType=bed&hgta_outFileName=SegmentalDups_GRCh38.bed

```bash
cut -f1-3 SegmentalDups_GRCh38.bed | sort -k1,1 -k2,2n | bedtools merge -i - | awk 'BEGIN {OFS="\t"} {print $0, "segmentaldup", "GRCh38"}' | awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/'> merged_SegmentalDups_GRCh38.bed
cut -f1-3 SegmentalDups_GRCh37.bed | sort -k1,1 -k2,2n | bedtools merge -i - | awk 'BEGIN {OFS="\t"} {print $0, "segmentaldup", "GRCh37"}' | awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/'> merged_SegmentalDups_GRCh37.bed
```

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


## Major Histocompatibility Complex (MHC) region
From :

GRCh37.p13
https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37.p13

GRCh38.p14
https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38.p14

### MHC region coordinates
chr6	28477797	33448354	MHC	GRCh37
chr6	28510120	33480577	MHC	GRCh38


## Telomeric and Centromeric regions dataset
Downloaded on 25/04/2025 from https://genome.ucsc.edu/cgi-bin/hgTables

ChromosomeBand_GRCh37.tsv from :
https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=2529613476_dkAUVDEoH74j8LaCc6nSM9DQngP5&clade=mammal&org=Human&db=hg19&hgta_group=map&hgta_track=cytoBand&hgta_table=0&hgta_regionType=genome&position=chr7%3A155%2C592%2C223-155%2C605%2C565&hgta_outputType=primaryTable&hgta_outFileName=ChromosomeBand_GRCh37.tsv


ChromosomeBand_GRCh38.tsv from :
https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=2529613476_dkAUVDEoH74j8LaCc6nSM9DQngP5&clade=mammal&org=Human&db=hg38&hgta_group=map&hgta_track=cytoBand&hgta_table=0&hgta_regionType=genome&position=chr7%3A155%2C592%2C223-155%2C605%2C565&hgta_outputType=primaryTable&hgta_outFileName=ChromosomeBand_GRCh38.tsv






### Formatting Code, example on GRCh37

#### Get first and last bands for each chromosome (after skipping header), then filtering out 'gvar' (acrocentric chromosome), keep only canonical chromosomes, and format the output with a "telomere" label.
```bash 
genome_version=GRCh37
awk 'NR > 1' ChromosomeBand_${genome_version}.tsv | sort -k1,1 -k2,2n | \
awk '
{
    chrom=$1
    if (chrom != prev_chrom) {
        if (NR > 2) print last_line
        print $0
        prev_chrom = chrom
    }
    last_line = $0
}
END {
    print last_line
}' | awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/' | grep -v gvar | cut -f1-3 | bedtools merge -i - | awk -v gv="${genome_version}" 'BEGIN {OFS="\t"} {print $0, "telomere", gv}' > telomere_${genome_version}.tsv
```
#### Extract centromeric regions (gieStain == "acen"), restrict to canonical chromosomes, and format with a "centromere" label.

```bash
grep acen ChromosomeBand_${genome_version}.tsv | awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/' | cut -f1-3 | bedtools merge -i - | awk -v gv="${genome_version}" 'BEGIN {OFS="\t"} {print $0, "centromere", gv}' > centromere_${genome_version}.tsv
```
#### Combine centromere and telomere regions into a unified, sorted file.

```bash
cat centromere_${genome_version}.tsv telomere_${genome_version}.tsv | sort -k1,1 -k2,2n > regions_${genome_version}.tsv
```
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


### Final Formatting

Problematic regions were then contatenated with previous annotations. In order to reduce redundancy, `bedtools merge` was performed on the final files. PAR1-2, XTR, MHC regions and the header was added back in manually.

For GRCh38
```bash
sort -k1,1 -k2,2n  Problematic_GRCh38_regions.tsv >  sorted_new_38.bed 
sort -k1,1 -k2,2n  regions_GRCh38.tsv > sorted_old_38.bed 

bedtools intersect -a sorted_old_38.bed -b sorted_new_38.bed  -v > unique_regions_38.bed 
cat sorted_old_38.bed  unique_regions_38.bed  | sort -k1,1 -k2,2n |  bedtools merge -i - | awk -F'\t' -v OFS='\t' '{$4="problematic_regions"; $5="GRCh38"; print}' > corrected_GRCh38.tsv
```

For GRCh37
```bash
sort -k1,1 -k2,2n  Problematic_GRCh37_regions.tsv >  sorted_new_37.bed 
sort -k1,1 -k2,2n  regions_GRCh37.tsv > sorted_old_37.bed 

bedtools intersect -a sorted_old_37.bed -b sorted_new_37.bed  -v > unique_regions_37.bed 
cat sorted_old_37.bed  unique_regions_37.bed  | sort -k1,1 -k2,2n |  bedtools merge -i - | awk -F'\t' -v OFS='\t' '{$4="problematic_regions"; $5="GRCh37"; print}' > corrected_GRCh37.tsv
```

To merge the two genome assemblies
```bash
cat corrected_GRCh38.tsv corrected_GRCh37.tsv >  Genome_Regions_data.tsv
```


