This is a nextflow module for annotating CNVs with Gene regions using VEP.

Requirements:
    - nextflow v24.X

    - python with polars installation; check the test.conf on how I configured in compute canada. For now I'm simply generating an environment inside the slurm job.

    - VEP container with a VEP cache

    - gnomad CNV files, the installation script is included as get_resources.sh

    - duckdb. The installation is as simple as:

        ```
        curl https://install.duckdb.org | sh
        ```

#TODO

 - Add process to ingest data from the merged PC and QS CNVs. The test data was generated from digCNV so we need something to handle the proper inputs
 - Hg19 compatibility?
