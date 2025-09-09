# ------------------------------------------------------------------------------
#
# Description:
# This script processes recurrent CNV data from an Excel file and exports
# CNV regions as BED files for both GRCh37 and GRCh38 genome versions.
#
# Main steps:
# 1. Load the recurrent CNV dataset from Excel.
# 2. Clean and standardize column names.
# 3. Parse genomic coordinates (Chr, Start, Stop) for GRCh37 and GRCh38.
# 4. Extract BED-format columns (Chr, Start, Stop).
# 5. Save BED files for each genome build.
# ------------------------------------------------------------------------------


library(readxl)
library(dplyr)
library(tidyr)

# Load Excel file
rCNV_df <- read_excel("~/flben/Git/CNV-Annotation/resources/rCNV/docs/recurrent_CNV_dataset.xlsx")

# Assign clean column names
colnames(rCNV_df) <- rCNV_df[1, ]

# Remove first row
rCNV_df <- rCNV_df[-1, ]

# Process each genome version separately
for(genome_version in c("GRCh37","GRCh38")){
  
  # Work on a copy, not the main dataframe
  tmp <- rCNV_df %>%
    tidyr::separate({{genome_version}}, into = c("Chr", "range"), sep = ":", remove = FALSE) %>%
    tidyr::separate(range, into = c("Start", "Stop"), sep = "-", remove = TRUE) %>%
    mutate(across(c(Start, Stop), as.numeric))
  
  # Select BED columns
  bed <- tmp[, c("Chr", "Start", "Stop")]
  
  # Output path
  out_file <- paste0("~/flben/recurrent_CNV/data/processed/clean/cnv_regions_", genome_version, ".bed")
  
  # Save BED file
  write.table(bed,
              file = out_file,
              sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}
