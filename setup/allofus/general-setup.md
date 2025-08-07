## General Set-Up Instructions


### Install DuckDB

The pipeline uses the CLI client of duckdb for creating parquet files and joining tables efficiently.  To install duckdb run:

```bash
curl https://install.duckdb.org | sh
export PATH="$HOME/.duckdb/cli/latest/:$PATH"
```

### Install Python dependencies 

```bash
pip install duckdb pandas matplotlib tqdm polars
```


### Install VEP

Our pipeline requires VEP 113 to annotate CNVs for transcripts. It can be installed from git. When prompted to cache files, type 'n' as we will do it ourselves in the next step.

```bash
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
git checkout release/113
perl INSTALL.pl
```



### Download Resources
The pipeline requires a VEP cache for annotations and auxiliary files from Gnomad.

For downloading the VEP cache for both Hg19 and 38 we use VEP version 113. To download the files use following commands:

```bash
curl -O https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz
tar xzf homo_sapiens_vep_113_GRCh38.tar.gz
curl -O https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh37.tar.gz
tar xzf homo_sapiens_vep_113_GRCh37.tar.gz
```

The tar commands will generate a directory called homo_sapiens. Gnomad auxiliary files should be stored there
```bash
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



### Modify Config

The project specific parameters and configurations can be set in the conf directory in a new config file. These configurations should be specific to the platform you plan on using. 

#### VEP
All processes requiring VEP are labeled with 'vep' in the process definitions.:

```yaml
     withLabel: vep {
        cpus = 16
        container = '/lustre06/project/6008022/All_user_common_folder/SOFTWARE/VEP/vep.sif'
        clusterOptions = '--mem-per-cpu 3500MB'
        time = '1h'
        }
```

#### Parameters
Paths to the relevent gnomad files and VEP cache must be provided:

```yaml
params {
    vep_cache      = "/lustre06/project/6008022/All_user_common_folder/SOFTWARE/VEP/cache
}

``` 