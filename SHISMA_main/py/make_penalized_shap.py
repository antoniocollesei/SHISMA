import os
import argparse
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.tree import DecisionTreeClassifier
from joblib import Parallel, delayed
import fasttreeshap
import importlib.metadata
from fast_borf import BorfBuilder
from fast_borf.pipeline.zero_columns_remover import ZeroColumnsRemover
from fast_borf.pipeline.reshaper import ReshapeTo2D
from fast_borf.pipeline.to_scipy import ToScipySparse
from fast_borf.xai.mapping import BagOfReceptiveFields
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
plt.style.use('default')
from constants import CUSTOM_CONFIG_A3, CUSTOM_CONFIG_A3_NO_DILATION, CUSTOM_CONFIG_A3_NO_DILATION_WINDOW_SIZE_2_3_4

# Argument parser
parser = argparse.ArgumentParser(description="Process time-series data with PPI filtering and BORF analysis.")
parser.add_argument("--time_series", required=True, help="Path to time-series data CSV file")
parser.add_argument("--ppi", required=True, help="Path to PPI data TSV file")
parser.add_argument("--cell_type", required=True, help="Cell type to analyze")
parser.add_argument("--output_folder", required=True, help="Output folder to store analysis results")
args = parser.parse_args()

# Load data
data = pd.read_csv(args.time_series, index_col=0).dropna(axis=0)
data.index = data.index.str.split('_').map(lambda x: (x[0], x[-1]))
data.index = pd.MultiIndex.from_tuples(data.index, names=["gene", "celltype"])
data_avg = data.groupby(level=["gene", "celltype"]).mean().reset_index()
data_avg["gene_celltype"] = data_avg["gene"] + "_" + data_avg["celltype"]
data_avg = data_avg.set_index("gene_celltype").drop(columns=["gene", "celltype"])

ppi_genes = pd.read_csv(args.ppi, sep="\t")
ppi_genes = pd.concat([ppi_genes["gene1"], ppi_genes["gene2"]]).unique()
data = data_avg.loc[data_avg.index.str.split('_').str[0].isin(ppi_genes)]

# Encoding
genes = np.array([name.split("_")[0] for name in list(data.index)])
cells = np.array([name.split("_")[1] for name in list(data.index)])
enc_genes, enc_cells = LabelEncoder(), LabelEncoder()
enc_genes.fit(genes)
enc_cells.fit(cells)
X, y_genes, y_cells = data.values[:, np.newaxis, :], enc_genes.transform(genes), enc_cells.transform(cells)

# BORF setup
builder = BorfBuilder(n_jobs=-2, configs=CUSTOM_CONFIG_A3_NO_DILATION_WINDOW_SIZE_2_3_4,
    pipeline_objects=[(ReshapeTo2D, dict(keep_unraveled_index=True)), (ZeroColumnsRemover, dict(axis=0)), (ToScipySparse, dict())])
borf = builder.build(X)
X_transformed = borf.fit_transform(X)

clf = DecisionTreeClassifier()
clf.fit(X_transformed, y_genes)
cell_type_indexes = np.where(enc_cells.inverse_transform(y_cells) == args.cell_type)[0]

shap_explainer = fasttreeshap.TreeExplainer(clf, algorithm='auto', n_jobs=-1)
shap_values = shap_explainer(X_transformed.toarray()[cell_type_indexes], check_additivity=False).values
average = np.mean(abs(shap_values), axis=0)
average_df = pd.DataFrame(average, columns=enc_genes.classes_.astype(str))

num_permutations, num_jobs = 100, -1
X_transformed_arr = X_transformed.toarray()[cell_type_indexes]

penalization_method = 'smash_to_zero_bonferroni'  # 'permutation' or 'maxT' or 'Westfall-Young'

# Permutation testing
if penalization_method == 'permutation':
  
  def compute_permutation(i):
    y_genes_perm = np.random.permutation(y_genes)
    model = DecisionTreeClassifier().fit(X_transformed, y_genes_perm)
    shap_explainer = fasttreeshap.TreeExplainer(model, algorithm='auto', n_jobs=num_jobs)
    shap_values_perm = shap_explainer(X_transformed_arr, check_additivity=False).values
    return np.einsum('ijk->jk', np.abs(shap_values_perm)) / shap_values_perm.shape[0]

  average_perm_list = np.array(Parallel(n_jobs=num_jobs)(delayed(compute_permutation)(i) for i in range(num_permutations)))
  p_values = pd.DataFrame(np.mean(average_perm_list >= average, axis=0))
  shap_values_pval_penalized = abs(average * (pd.DataFrame(1-p_values)))
  shap_values_pval_penalized.columns = enc_genes.classes_.astype(str)

elif penalization_method == 'smash_to_zero_bonferroni':
  
  penalization_threshold = 0.05
  num_features = average.shape[1]  # Number of features (genes)
  bonferroni_threshold = penalization_threshold / num_features  # Adjusted threshold

  def compute_permutation(i):
    y_genes_perm = np.random.permutation(y_genes)
    model = DecisionTreeClassifier().fit(X_transformed, y_genes_perm)
    shap_explainer = fasttreeshap.TreeExplainer(model, algorithm='auto', n_jobs=num_jobs)
    shap_values_perm = shap_explainer(X_transformed_arr, check_additivity=False).values
    return np.einsum('ijk->jk', np.abs(shap_values_perm)) / shap_values_perm.shape[0]

  average_perm_list = np.array(Parallel(n_jobs=num_jobs)(delayed(compute_permutation)(i) for i in range(num_permutations)))
  p_values = pd.DataFrame(np.mean(average_perm_list >= average, axis=0))
  penalization_mask = (p_values <= bonferroni_threshold).astype(int)
  shap_values_pval_penalized = abs(average * penalization_mask)
  shap_values_pval_penalized.columns = enc_genes.classes_.astype(str)

elif penalization_method == 'smash_to_zero_fdr':
  
  penalization_threshold = 0.05  # Original significance level
  num_features = average.shape[1]  # Number of features (genes)

  def compute_permutation(i):
    y_genes_perm = np.random.permutation(y_genes)
    model = DecisionTreeClassifier().fit(X_transformed, y_genes_perm)
    shap_explainer = fasttreeshap.TreeExplainer(model, algorithm='auto', n_jobs=num_jobs)
    shap_values_perm = shap_explainer(X_transformed_arr, check_additivity=False).values
    return np.einsum('ijk->jk', np.abs(shap_values_perm)) / shap_values_perm.shape[0]

  average_perm_list = np.array(Parallel(n_jobs=num_jobs)(delayed(compute_permutation)(i) for i in range(num_permutations)))
  p_values = pd.DataFrame(np.mean(average_perm_list >= average, axis=0))

  p_values_corrected = multipletests(p_values.values.flatten(), alpha=penalization_threshold, method='fdr_bh')[1]
  p_values_corrected = p_values_corrected.reshape(p_values.shape)  # Reshape to match original format

  # Apply FDR correction threshold
  penalization_mask = (p_values_corrected <= penalization_threshold).astype(int)

  # Penalize SHAP values
  shap_values_pval_penalized = pd.DataFrame(abs(average * penalization_mask), columns=enc_genes.classes_.astype(str))

os.makedirs(args.output_folder, exist_ok=True)
shap_values_pval_penalized.to_csv(f"{args.output_folder}/shap_values_{args.cell_type}_pvalpenalized.csv", index=False)

### Create input files for hierarchical-hotnet
shap_values_pval_penalized = shap_values_pval_penalized.abs()
genes = shap_values_pval_penalized.columns

for i in range(len(shap_values_pval_penalized)):
  row = shap_values_pval_penalized.iloc[i]
  if row.sum() != 0:
    row.T.to_csv(f"{args.output_folder}/shapelet_data/scores_{i}.tsv", header=False, sep="\t")

# Create a gene-index file
gene2index = pd.Series(index=genes, data=range(len(genes)))

# Load edges
edges = pd.read_csv(args.ppi, sep="\t")
network_name = os.path.splitext(os.path.basename(args.ppi))[0]

# Keep edges that are in genes
edges = edges[edges['gene1'].isin(genes) & edges['gene2'].isin(genes)]

# Map genes to indexes
edges['gene1'] = edges['gene1'].map(gene2index)
edges['gene2'] = edges['gene2'].map(gene2index)

# Print to csv only gene1, gene2
edge_list = edges[['gene1', 'gene2']]
edge_list.to_csv(f"{args.output_folder}/{network_name}_edge_list.tsv", header=False, sep="\t", index=False)

index2gene = pd.Series(data=genes, index=range(len(genes)))
index2gene.to_csv(f"{args.output_folder}/{network_name}_index_gene.tsv", header=False, sep="\t")

# Export mappings
mapper = BagOfReceptiveFields(borf)
mapper.build(X)

# Create a PDF file
with PdfPages(f"{args.output_folder}/mappings.pdf") as pdf:
  for i in range(average.shape[0]):
    plt.figure()
    plt.plot(mapper[i].word_array, marker="o")
    plt.title(f'Shapelet number {i}')
    pdf.savefig()  # Save the current figure into the PDF
    plt.close()  # Close the figure to free memory

print("Processing complete.")
