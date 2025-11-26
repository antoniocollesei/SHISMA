# needed:
# - time labeled seurat objects
# - patient_id metadata column in each seurat object
# - cell_type metadata column in each seurat object

# load libraries
library(tidyverse)
library(Seurat)

# create folder with output files, if not exists
dir.create("data_ready", showWarnings = FALSE)

d0 <- readRDS("d0.rds")
d1 <- readRDS("d1.rds")
d2 <- readRDS("d2.rds")
d5 <- readRDS("d5.rds")
d9 <- readRDS("d9.rds")
d15 <- readRDS("d15.rds")

list_seurat_obj <- ls()
list_seurat_obj <- list_seurat_obj[order(as.numeric(gsub("d", "", list_seurat_obj)))] # sort list_seurat_obj by timepoint

# select assay of reference
assay_selected <- "SCT"

# loop through each seurat object, create dataframe with mean expression per gene per patient per celltype, and merge all dataframes to create the time series dataframe
for (seurat_obj in list_seurat_obj) {
  tmp_obj <- get(seurat_obj)

  # create dataframe from tmp_obj
  tmp_df <- as.data.frame(tmp_obj@assays[[assay_selected]]$data)

  colnames(tmp_df) <- paste0(
    tmp_obj$patient_id,
    "_",
    tmp_obj$cell_type
  ) %>% 
    make.names()

  # group by patient and celltype, and do the mean
  tmp_df_mean <- t(apply(tmp_df, 1, function(x) tapply(x, colnames(tmp_df), mean))) %>% as.data.frame()

  # group by identical colnames and mean expression
  tmp_df_mean <- rownames_to_column(tmp_df_mean, var = "gene")

  tmp_df_long <- tmp_df_mean %>%
    pivot_longer(cols = -gene, names_to = "colname", values_to = "value")

  # paste together gene and colname (= patient_id + celltype) into new column
  tmp_df_long$gene_patient_celltype <- paste(
    tmp_df_long$gene,
    tmp_df_long$colname,
    sep = "_"
  )

  # remove gene and colname
  tmp_df_long <- tmp_df_long %>%
    select(gene_patient_celltype, value)

  # append to dataframe by joining with gene_patient_celltype
  if (seurat_obj == list_seurat_obj[1]) {
    final_df <- tmp_df_long
  } else {
    final_df <- full_join(final_df, tmp_df_long, by = "gene_patient_celltype")
  }

  # print progress
  print(paste("Finished", seurat_obj))
}

colnames(final_df) <- c("gene_patient_celltype", list_seurat_obj)

# replace NA with 0
final_df[is.na(final_df)] <- 0

# write final_df to csv
write.csv(final_df, "data_ready/gene_patient_celltype_normalized.csv", row.names = FALSE)
