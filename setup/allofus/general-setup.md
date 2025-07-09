## General Set-Up Instructions


### Install DuckDB

The pipeline uses the CLI client of duckdb for creating parquet files and joining tables efficiently.  To install duckdb run:

```
curl https://install.duckdb.org | sh
```

### Download Resources
The pipeline requires a VEP cache for annotations and auxiliary files from Gnomad.

For downloading the VEP cache for both Hg19 and 38 we use VEP version 113. To download the files use following commands:

```
curl -O https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz
tar xzf homo_sapiens_vep_113_GRCh38.tar.gz
curl -O https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh37.tar.gz
tar xzf homo_sapiens_vep_113_GRCh37.tar.gz
```

The tar commands will generate a directory called homo_sapiens. Gnomad auxiliary files can be stored there or elsewhere using the following script and replacing the "dir" variable:
```
#Getting allele frequency scores for Structural Variants (SVs)
#https://gnomad.broadinstitute.org/news/2023-11-v4-copy-number-variants/

dir="gnomad_resources"

#V4 Hg38
curl "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/exome_cnv/gnomad.v4.1.cnv.all.vcf.gz"  -o "$dir/gnomad.v4.1.cnv.all.vcf.gz" &&

                                                                        gunzip -c  "$dir/gnomad.v4.1.cnv.all.vcf.gz" |
                                                                        bgzip > "$dir/gnomad.v4.1.cnv.all.vcf.bgz"   &&
                                                                        tabix -p vcf "$dir/gnomad.v4.1.cnv.all.vcf.bgz"

curl "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz"  -o  "$dir/gnomad.v4.1.sv.sites.vcf.gz" &&
                                                                        gunzip -c  "$dir/gnomad.v4.1.sv.sites.vcf.gz" |
                                                                        bgzip > "$dir/gnomad.v4.1.sv.sites.vcf.bgz"   &&
                                                                        tabix -p vcf "$dir/gnomad.v4.1.sv.sites.vcf.bgz"

#V2 Hg19
curl "https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz" -o "$dir/gnomad_v2.1_sv.sites.vcf.gz" &&
                                                                        gunzip -c "$dir/gnomad_v2.1_sv.sites.vcf.gz" |
                                                                        bgzip > "$dir/gnomad.v2.1.sv.sites.vcf.bgz"   &&
                                                                        tabix -p vcf "$dir/gnomad.v2.1.sv.sites.vcf.bgz"
#Getting genetic Constraint metrics (pLi and LOEUF)
#https://gnomad.broadinstitute.org/help/constraint

curl "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv" -o "$dir/gnomad.v4.1.constraint_metrics.tsv"

```




### Build VEP Container
The pipeline requires version 113 of VEP. Depending on the platform different containerization software can be used to provide the VEP executable. Regardless, Nextflow supports both docker and apptainer for mounting an image. First, it is recommended to use the following container from dockerhub:

```
 ensemblorg/ensembl-vep:release_113.4
```

For use on compute canada, compute nodes cannot access the internet, therefore we use a local image built using apptainer.

```
apptainer pull vep.sif docker://ensemblorg/ensembl-vep:release_113.4
```




### Modify Config

The project specific parameters and configurations can be set in the conf directory in a new config file. These configurations should be specific to the platform you plan on using. 

#### VEP
All processes requiring VEP are labeled with 'vep' in the process definitions. In the process field the absolute path to the container should be provided:

```
     withLabel: vep {
        cpus = 16
        container = '/lustre06/project/6008022/All_user_common_folder/SOFTWARE/VEP/vep.sif'
        clusterOptions = '--mem-per-cpu 3500MB'
        time = '1h'
        }
```

#### Parameters
Paths to the relevent gnomad files and VEP cache must be provided:

```
params {
    vep_cache      = "/lustre06/project/6008022/All_user_common_folder/SOFTWARE/VEP/cache"
    genome_regions = "${projectDir}/resources/Genome_Regions/Genome_Regions_data.tsv"
    
    //HG38 specific gnomad files
    gnomad_AF          = "/lustre06/project/6008022/All_user_common_folder/SOFTWARE/VEP/cache/ressources_gnomAD/gnomad.v4.1.sv.sites.vcf.bgz"
    gnomad_constraints = "/lustre06/project/6008022/All_user_common_folder/SOFTWARE/VEP/cache/ressources_LOEUF/gnomad.v4.1.constraint_metrics.tsv"
}

``` 