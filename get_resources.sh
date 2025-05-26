#!/bin/bash
module load htslib/1.19

#https://gnomad.broadinstitute.org/news/2023-11-v4-copy-number-variants/

curl "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/exome_cnv/gnomad.v4.1.cnv.all.vcf.gz" | 
                                                                                                     gunzip -c | 
                                                                        bgzip > gnomad.v4.1.cnv.all.vcf.bgz   &&
                                                                        tabix -p vcf gnomad.v4.1.cnv.all.vcf.bgz 

curl "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz" |
                                                                                                     gunzip -c  | 
                                                                        bgzip > gnomad.v4.1.sv.sites.vcf.bgz   &&
                                                                        tabix -p vcf gnomad.v4.1.sv.sites.vcf.bgz


curl "https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz" |                                                                                                      
											      gunzip -c  | 
                                                                        bgzip > gnomad.v2.1.sv.sites.vcf.bgz   &&
                                                                        tabix -p vcf gnomad.v2.1.sv.sites.vcf.bgz
