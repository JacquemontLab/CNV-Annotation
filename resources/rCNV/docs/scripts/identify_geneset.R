# ------------------------------------------------------------------------------
# 
# Description:
# This script processes recurrent CNV data from an Excel file and generates a
# unified geneset per rCNV entry for both GRCh37 and GRCh38 genome versions.
#
# Main steps:
# 1. Load the recurrent CNV dataset from Excel.
# 2. Clean and standardize column names.
# 3. Parse genomic coordinates (Chr, Start, Stop) for GRCh37 and GRCh38.
# 4. Create unique recurrent CNV identifiers per genome build.
# 5. Load overlapping annotated genes for each genome version.
# 6. Merge CNVs with gene annotations and collapse into non-redundant genesets.
# 7. Attach the geneset for each CNV (for both GRCh37 and GRCh38).
# 8. Filter CNVs to retain only those with genesets in both genome versions.
# 9. Save the final table (`geneset_per_rCNV.tsv`) for downstream analyses.
#
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

# Work on a copy, not the main dataframe
# Process each genome version separately
for(genome_version in c("GRCh37","GRCh38")){
  tmp <- rCNV_df %>%
    tidyr::separate({{genome_version}}, into = c("Chr", "range"), sep = ":", remove = FALSE) %>%
    tidyr::separate(range, into = c("Start", "Stop"), sep = "-", remove = TRUE) %>%
    mutate(across(c(Start, Stop), as.numeric))
  
  # dynamically assign column name
  tmp[[paste0("recurrent_ID_", genome_version)]] <- paste(tmp$Chr, tmp$Start, tmp$Stop, sep="_")
  
  
  
  overlapping_gene <- paste0("~/flben/recurrent_CNV/data/processed/clean/dummy_annotated_filtered_genes_", genome_version, ".tsv")
  
  gene_set = read.delim(overlapping_gene,header=T)
  gene_set[[paste0("recurrent_ID_", genome_version)]] = paste(gene_set$Chr,gene_set$Start,gene_set$End,sep="_")
  
  gene_set = gene_set[!duplicated(paste(gene_set[[paste0("recurrent_ID_", genome_version)]],gene_set$Gene,sep="_")),]
  
  
  
  merge_df = merge.data.frame(x = tmp, y = gene_set, by=paste0("recurrent_ID_", genome_version))
  
  collapsed <- merge_df %>%
    group_by(rCNV_ID) %>%
    summarise(Genes = paste(unique(Gene), collapse = ","), .groups = "drop")
  
  rownames(collapsed) = collapsed$rCNV_ID
  
  rCNV_df[[paste0("geneset_", genome_version)]] = collapsed[rCNV_df$rCNV_ID, ]$Genes
}

rCNV_df = rCNV_df[!is.na(rCNV_df$geneset_GRCh38) & !is.na(rCNV_df$geneset_GRCh37),]

write.table(rCNV_df,
            file = "~/flben/recurrent_CNV/data/processed/geneset_per_rCNV.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

