library(readxl)
library(dplyr)
library(tidyr)

# Load Excel file
rCNV_df <- read_excel("~/flben/Git/CNV-DB-Builder/resources/rCNV/docs/recurrent_CNV_dataset.xlsx")

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
  
  tmp$recurrent_ID <- paste(tmp$Chr, tmp$Start, tmp$Stop, sep="_")
  
  # Select BED columns
  bed <- tmp[, c("Chr", "Start", "Stop")]
  
  # Output path
  out_file <- paste0("~/flben/recurrent_CNV/data/processed/clean/cnv_regions_", genome_version, ".bed")
  
  # Save BED file
  write.table(bed,
              file = out_file,
              sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}
