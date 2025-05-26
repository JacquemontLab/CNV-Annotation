#!/bin/bash
#SBATCH --job-name=nextflow_job_test    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G                     # Job memory request
#SBATCH --time=01:00:00               # Time limit hrs:min:sec
#SBATCH --output=nextflow_test_%j.log   # Standard output and error log
#SBATCH --account=rrg-jacquese        #group account


module load java

NXF_VER=25.04.2
~/nextflow run main.nf -c conf/ccdb.config -with-report
