#!/bin/bash
# First argument is the .tsv file to annotate

#SBATCH --job-name=cnvDB_Buider    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=3500MB                     # Job memory request
#SBATCH --time=01:00:00               # Time limit hrs:min:sec
#SBATCH --output=cnvDB_Buider_%j.log   # Standard output and error log
#SBATCH --account=rrg-jacquese        #group account


module load java
module load python/3.13.2
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index polars


NXF_VER=25.04.2
export NXF_OFFLINE=true
/lustre06/project/6008022/All_user_common_folder/SOFTWARE/Nextflow/nextflow-25.04.2-dist run main.nf --cnvs $1 \
    --genome_version "GRCh37" \
    --cohort_tag "cnvDB_Buider" \
    -c setup/ccdb/ccdb.config \
    -with-report results/SPARK_microarray/docs/report.html \
    -resume
