
library(LEAP)
library(dplyr)

data <- read.csv("../../SHISMA_main/temporal_data_with_patient_ready_normalized_full_genes.csv")
data_safe <- data

# keep rows that contain Bcells in rownames
data <- data[grep("Bcells", data$gene_patient_celltype), ]
data %>% head()

# rownames containing the same string before the first underscore get averaged
data$gene_patient_celltype <- sapply(strsplit(data$gene_patient_celltype, "_"), `[`, 1)
data_agg <- aggregate(. ~ gene_patient_celltype, data = data, FUN = mean)
data_agg %>% head()

# load ppi for filtering
ppi <- read.csv("../../SHISMA_main/string_is_0.7_ev_reactome.tsv", sep = "\t")
data_agg <- data_agg[data_agg$gene_patient_celltype %in% unique(c(ppi$gene1, ppi$gene2)), ]

rownames(data_agg) <- data_agg$gene_patient_celltype
data_agg <- data_agg[, -1]
data <- data_agg

data %>% head() 

# launch LEAP
output <- MAC_counter(data=data, max_lag_prop=1/3, file_name="real_data_Bcells", symmetric=TRUE)
output <- output %>% as.data.frame()
colnames(output) <- c("Correlation", "Lag", "Gene1", "Gene2")
output$Gene1 <- rownames(data)[output$Gene1]
output$Gene2 <- rownames(data)[output$Gene2]

rownames(output) <- NULL
output %>% head()

write.csv(output, file="real_data_Bcells_LEAP_output.csv", row.names=FALSE)
