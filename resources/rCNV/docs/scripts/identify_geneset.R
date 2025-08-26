library(readxl)
library(dplyr)
library(tidyr)

# Load Excel file
rCNV_df <- read_excel("~/flben/Git/CNV-DB-Builder/resources/rCNV/docs/recurrent_CNV_dataset.xlsx")

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

