library(tidyverse)

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
  tmp_df <- as.data.frame(tmp_obj@assays$RNA$scale.data)
  colnames(tmp_df) <- tmp_obj$cell_type

  # group by identical colnames and mean expression
  tmp_df$rowname <- rownames(tmp_df)

  tmp_df_long <- tmp_df %>%
    pivot_longer(cols = -rowname, names_to = "colname", values_to = "value")

  tmp_df_mean <- tmp_df_long %>%
    group_by(rowname, colname) %>%
    summarise(mean_value = mean(value), .groups = "drop")

  # paste together rowname and colname into new column
  tmp_df_mean$gene_cell_type <- paste(
    tmp_df_mean$rowname,
    tmp_df_mean$colname,
    sep = "_"
  )

  # remove spaces from gene_cell_type
  tmp_df_mean$gene_cell_type <- gsub(" ", "_", tmp_df_mean$gene_cell_type)

  # remove rowname and colname
  tmp_df_mean <- tmp_df_mean %>%
    select(gene_cell_type, mean_value)

  # append to dataframe by joining with gene_cell_type
  if (seurat_obj == "d0_obj") {
    final_df <- tmp_df_mean
  } else {
    final_df <- full_join(final_df, tmp_df_mean, by = "gene_cell_type")
  }

  # print progress
  print(paste("Finished", seurat_obj))
}

# get only what's before the first underscore of list_seurat_obj
list_seurat_obj <- gsub("_.*", "", list_seurat_obj)

colnames(final_df) <- c("gene_cell_type", list_seurat_obj)

final_df <- final_df %>%
  select(gene_cell_type, d0, d1, d2, d5, d9, d15)

# remove tmp objects
rm(tmp_df, tmp_df_long, tmp_df_mean, tmp_obj)

# write final_df to csv
write.csv(final_df, "temporal_data_ready_normalized.csv", row.names = FALSE)
