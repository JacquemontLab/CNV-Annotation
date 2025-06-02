#!/bin/bash
#SBATCH --job-name=SPARK_XHMM_DB    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G                     # Job memory request
#SBATCH --time=01:00:00               # Time limit hrs:min:sec
#SBATCH --output=SPARK_XHMM_DB_%j.log   # Standard output and error log
#SBATCH --account=rrg-jacquese        #group account


module load java


module load python/3.13.2
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index polars


NXF_VER=25.04.2
~/nextflow run main.nf  --cnvs test-data/XHMM_allCNV_SPARK_order.tsv \
--genome_version "GRCh38" \
--cohort_tag "SPARK_XHMM" \
-c conf/ccdb.config \
-with-report

