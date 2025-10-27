import os, json, glob
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.preprocessing import LabelEncoder, StandardScaler, FunctionTransformer

def aggregate_clusters(results_dir, output_file):
    aggregated_data = []
    min_subnetwork_size = 3  # Exclude single-gene subnetworks and couples

    # Find all cluster result files
    cluster_files = glob.glob(os.path.join(results_dir, "clusters_*.tsv"))
    
    for file_path in cluster_files:
        with open(file_path, "r") as f:
            lines = f.readlines()
        
        # Extract network and score number from filename
        filename = os.path.basename(file_path)
        parts = filename.replace("clusters_", "").replace(".tsv", "").split("_scores_")
        if len(parts) != 2:
            continue
        network = parts[0]
        score_number = parts[1]  # Capture score number
        
        # Find clusters section and p-value
        cluster_start = next((i for i, line in enumerate(lines) if line.strip() == "# Clusters:"), None)
        p_value = None
        for line in lines:
            if line.startswith("# p-value:"):
                p_value = float(line.strip().split(": ")[1])  # Convert to float for sorting
                break
        
        if cluster_start is None:
            continue
        
        # Read clusters, excluding single-gene subnetworks
        for line in lines[cluster_start + 1:]:
            genes = line.strip().split("\t")
            if len(genes) >= min_subnetwork_size:
                aggregated_data.append([network, score_number, p_value, ",".join(genes)])
    
    # Create DataFrame
    df = pd.DataFrame(aggregated_data, columns=["Network", "Shapelet Number", "P-Value", "Genes"])
    
    # Sort by p-value (ascending order)
    df = df.sort_values(by="P-Value")
    df.reset_index(drop=True, inplace=True)
    
    # Save aggregated results
    df.to_csv(output_file, sep="\t", index=True)
    print(f"Aggregated results saved to {output_file}")
    
    return df

# Argument parser
parser = argparse.ArgumentParser(description="Process time-series data with PPI filtering and BORF analysis.")
parser.add_argument("--time_series", required=True, help="Path to time-series data CSV file")
parser.add_argument("--cell_type", required=True, help="Chosen cell-type")
parser.add_argument("--results_dir", required=True, help="Path to hierarchical hotnet results directory")
parser.add_argument("--output_folder", required=True, help="Path to output directory")
args = parser.parse_args()

# Load data
data = pd.read_csv(args.time_series, index_col=0).dropna(axis=0)
data.index = data.index.str.split('_').map(lambda x: (x[0], x[-1]))
data.index = pd.MultiIndex.from_tuples(data.index, names=["gene", "celltype"])
data_avg = data.groupby(level=["gene", "celltype"]).mean().reset_index()
data_avg["gene_celltype"] = data_avg["gene"] + "_" + data_avg["celltype"]
data_avg = data_avg.set_index("gene_celltype").drop(columns=["gene", "celltype"])
data = data_avg

# Encoding
genes = np.array([name.split("_")[0] for name in list(data.index)])
cells = np.array([name.split("_")[1] for name in list(data.index)])
enc_genes, enc_cells = LabelEncoder(), LabelEncoder()
enc_genes.fit(genes)
enc_cells.fit(cells)
X, y_genes, y_cells = data.values[:, np.newaxis, :], enc_genes.transform(genes), enc_cells.transform(cells)

# Results
results_directory = args.results_dir
output_file = f"{args.output_folder}/aggregated_clusters_{args.output_folder}.tsv"

# Ensure the output directory exists
os.makedirs(args.output_folder, exist_ok=True)
print(f"Aggregating clusters from {results_directory} into {output_file}...")
results = aggregate_clusters(results_directory, output_file)

output_plot_dir = f"{args.output_folder}/plots"
os.makedirs(output_plot_dir, exist_ok=True)
print(f"Printing results into {output_plot_dir}...")

# Plotting only significant clusters
results["P-Value"] = results["P-Value"].astype(float)
results_filt = results[results["P-Value"] < 0.05]

# For each row of results_filt, find rows of data containing the genes in the cluster and the chosen cell type
for index, row in results_filt.iterrows():
  genes = row["Genes"].split(",")
  genes = [gene.strip() for gene in genes]
  cell_type = args.cell_type
  data_filtered = data.loc[data.index.str.split('_').str[0].isin(genes) & (data.index.str.split('_').str[1] == cell_type)]
  
  # Calculate mean and standard error for each gene
  data_filtered_mean = data_filtered.groupby(data_filtered.index.str.split('_').str[0]).mean()
  data_filtered_sem = data_filtered.groupby(data_filtered.index.str.split('_').str[0]).sem()
  
  # Plotting
  plt.figure(figsize=(10, 6))
  for gene in data_filtered_mean.index:
    plt.plot(data_filtered_mean.columns, data_filtered_mean.loc[gene], label=gene)
    plt.fill_between(
      data_filtered_mean.columns,
      data_filtered_mean.loc[gene] - data_filtered_sem.loc[gene],
      data_filtered_mean.loc[gene] + data_filtered_sem.loc[gene],
      alpha=0.2
    )
  
  plt.title(f"Cluster: {row['Network']}")
  plt.xlabel("Samples")
  plt.ylabel("Mean Expression (within cell type)")
  plt.legend(title="Genes", bbox_to_anchor=(1.05, 1), loc='upper left')
  plt.tight_layout()
  
  # Save the plot
  plot_file = os.path.join(output_plot_dir, f"cluster_{index}.pdf")
  plt.savefig(plot_file)
  plt.close()