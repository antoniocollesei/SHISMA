# R version: 4.4.2

library(tidyverse)

# create folder with output files, if not exists
dir.create("data_ready", showWarnings = FALSE)

d0_obj <- readRDS("rds/d0.rds")
d1_obj <- readRDS("rds/d1.rds")
d2_obj <- readRDS("rds/d2.rds")
d5_obj <- readRDS("rds/d5.rds")
d9_obj <- readRDS("rds/d9.rds")
d15_obj <- readRDS("rds/d15.rds")

list_seurat_obj <- ls()

for (seurat_obj in list_seurat_obj) {
  tmp_obj <- get(seurat_obj)

  # create dataframe from tmp_obj
  tmp_df <- as.data.frame(tmp_obj@assays$RNA$data)

  colnames(tmp_df) <- paste0(
    substr(colnames(tmp_obj), 1, 2),
    "_",
    tmp_obj$cell_type
  )

  tmp_df_mean <- t(apply(tmp_df, 1, function(x) tapply(x, colnames(tmp_df), mean))) %>% as.data.frame()

  # group by identical colnames and mean expression
  tmp_df_mean <- rownames_to_column(tmp_df_mean, var = "gene")

  tmp_df_long <- tmp_df_mean %>%
    pivot_longer(cols = -gene, names_to = "colname", values_to = "value")

  # paste together gene and colname into new column
  tmp_df_long$gene_patient_celltype <- paste(
    tmp_df_long$gene,
    tmp_df_long$colname,
    sep = "_"
  ) %>% gsub(" ", "", .)

  # remove gene and colname
  tmp_df_long <- tmp_df_long %>%
    select(gene_patient_celltype, value)

  # append to dataframe by joining with gene_patient_celltype
  if (seurat_obj == "d0_obj") {
    final_df <- tmp_df_long
  } else {
    final_df <- full_join(final_df, tmp_df_long, by = "gene_patient_celltype")
  }

  # print progress
  print(paste("Finished", seurat_obj))
}

# get only what's before the first underscore of list_seurat_obj
list_seurat_obj <- gsub("_.*", "", list_seurat_obj)

colnames(final_df) <- c("gene_patient_celltype", list_seurat_obj)

final_df <- final_df %>%
  select(gene_patient_celltype, d0, d1, d2, d5, d9, d15)

final_df[is.na(final_df)] <- 0

# keep everything before underscore
gsub("_.*", "", final_df$gene_patient_celltype) %>% unique()

# remove tmp objects
rm(tmp_df, tmp_df_long, tmp_df_mean, tmp_obj)

# write final_df to csv
write.csv(final_df, "temporal_data_with_patient_ready_normalized_full_genes.csv", row.names = FALSE)
