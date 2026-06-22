import numpy as np
import pandas as pd
import shap, fasttreeshap
from sklearn.base import clone
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.linear_model import RidgeClassifier, LogisticRegression
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from lightgbm import LGBMClassifier
from sklearn.preprocessing import LabelEncoder
from fast_borf import BorfBuilder
from fast_borf.pipeline.zero_columns_remover import ZeroColumnsRemover
from fast_borf.pipeline.reshaper import ReshapeTo2D
from fast_borf.pipeline.to_scipy import ToScipySparse
from fast_borf.xai.mapping import BagOfReceptiveFields
import matplotlib.pyplot as plt
from joblib import Parallel, delayed

import sys
from constants import CUSTOM_CONFIG_A3_NO_DILATION_WINDOW_SIZE_3

import time
import itertools
from scipy.stats import kendalltau
from joblib import Parallel, delayed

import warnings
warnings.filterwarnings('ignore')

import logging, os
logging.getLogger().setLevel(logging.ERROR)
plt.style.use('default')


num_permutations = 10
num_cores = -1
penalization_method = "bonferroni"
penalization_threshold = 0.05
output_folder = "./results"


def generate_shisma_data(n_genes=200, n_planted=20, n_timepoints=6, n_pats=5, n_celltypes=10):
    np.random.seed(42)
    
    cell_types = np.array([f"CT{i}" for i in range(n_celltypes)])
    target_ct = "CT0"
    padding = len(str(n_genes - 1))
    genes = np.array([f"G{i:0{padding}d}" for i in range(n_genes)])
    planted_genes = genes[:n_planted]
    
    # ==========================================
    # 1. PPI GENERATION
    # ==========================================
    edges = [[planted_genes[i], planted_genes[j]] for i in range(n_planted) for j in range(i+1, n_planted)]
    for _ in range(n_genes * 2): 
        g1, g2 = np.random.choice(genes, 2)
        if g1 != g2: edges.append([g1, g2])
    ppi_df = pd.DataFrame(edges, columns=["gene1", "gene2"]).drop_duplicates()

    # ==========================================
    # 2. REALISTIC BACKGROUND EXPRESSION (Negative Binomial)
    # ==========================================
    print("Generating realistic Negative Binomial gene expression...")
    n_samples = n_genes * n_celltypes * n_pats
    
    # Create aligned metadata grids for (Gene, CellType, pat)
    mesh_genes, mesh_cts, mesh_pats = np.meshgrid(np.arange(n_genes), np.arange(n_celltypes), np.arange(n_pats), indexing='ij')
    flat_gene_idx = mesh_genes.flatten()
    flat_ct_idx = mesh_cts.flatten()
    flat_pats = mesh_pats.flatten()
    
    # A) Gene-Specific Baselines: 
    # Some genes are naturally highly expressed, others are nearly silent (Gamma distribution)
    gene_baselines = np.random.gamma(shape=2.0, scale=1.5, size=n_genes)
    
    # B) Cell-Type Specific Variance: 
    # Cell types express the same genes at slightly different levels (Log-Normal multipliers)
    ct_multipliers = np.random.lognormal(mean=0.0, sigma=0.5, size=(n_genes, n_celltypes))
    
    # Calculate the base expected mean (lambda) for each specific row
    base_means = gene_baselines[flat_gene_idx] * ct_multipliers[flat_gene_idx, flat_ct_idx]
    
    # Expand across timepoints and add time-step biological noise
    mu_matrix = np.tile(base_means[:, np.newaxis], (1, n_timepoints))
    biological_noise = np.random.lognormal(mean=0.0, sigma=0.1, size=(n_samples, n_timepoints))
    mu_matrix = mu_matrix * biological_noise
    
    # C) Sample Counts: 
    # Drawing from a Poisson with Gamma-distributed means = Negative Binomial distribution
    raw_counts = np.random.poisson(lam=mu_matrix)
    
    # D) Zero-Inflation (Dropouts):
    # Simulate technical dropouts (lower expressed genes are more likely to drop to exactly 0)
    dropout_probs = np.exp(-0.1 * mu_matrix)
    dropouts = np.random.binomial(n=1, p=dropout_probs)
    raw_counts = raw_counts * (1 - dropouts)
    
    # E) Standard Transcriptomics Normalization (Log1p):
    # Converts integer counts to stable, continuous values (standard for ML pipelines)
    normalized_expression = np.log1p(raw_counts)
    
    # ==========================================
    # 3. OPTIONAL: PLANTING YOUR SIGNAL
    # ==========================================
    # If you ever want to inject your signal back into this realistic baseline, 
    # you would do it here on the continuous `normalized_expression` matrix:
    #
    is_target_ct = (cell_types[flat_ct_idx] == target_ct)
    is_planted_gene = np.isin(genes[flat_gene_idx], planted_genes)
    plant_mask = is_target_ct & is_planted_gene
    normalized_expression[plant_mask, 3] += 2.0
    normalized_expression[plant_mask, 4] -= 1.5 
    normalized_expression[plant_mask, 5] += 2.0
    
    # ==========================================
    # 4. DATAFRAME FORMATTING
    # ==========================================
    df_synth = pd.DataFrame(normalized_expression, columns=[f"T{t}" for t in range(n_timepoints)])
    
    flat_genes = genes[flat_gene_idx]
    flat_cts = cell_types[flat_ct_idx]
    
    df_synth.index = [f"{g}_pat{r}_{ct}" for g, ct, r in zip(flat_genes, flat_cts, flat_pats)]
    
    return df_synth, ppi_df, list(planted_genes), target_ct


df_synth, ppi_synth, target_genes, target_ct = generate_shisma_data(
    n_genes=300, 
    n_planted=30, 
    n_timepoints=6, 
    n_pats=6,        
    n_celltypes=5

df_synth.head(), df_synth.shape

# ==========================================
# 1. DATA PREP & BORF TRANSFORMATION
# ==========================================
print("1. Filtering for PPI genes & Handling Replicates...")

# Extract unique genes present anywhere in the PPI network
ppi_genes = np.unique(ppi_synth[['gene1', 'gene2']].values.astype(str))

# Parse the index into a MultiIndex (safely handling underscores in gene names)
# Assuming format: Gene_Celltype_Replicate OR Gene_Patient_Celltype
# Adjust the names below to match your actual string structure
index_dims = df_synth.index.str.rsplit('_', n=2, expand=True)
index_dims.names = ["gene", "patient", "celltype"] 
df_synth.index = index_dims

# 1a. Early PPI Filtering
df_synth = df_synth[df_synth.index.get_level_values("gene").isin(ppi_genes)]
print(f"  -> Reduced to {len(df_synth.index.get_level_values('gene').unique())} valid PPI genes.")

# 1b. Handle Replicates/Patients (Toggle Option A or B)
# ---------------------------------------------------------
# OPTION A: Average across replicates/patients (Collapses the 3rd level)
df_processed = df_synth#.groupby(level=["gene", "celltype"]).mean()
print("  -> Averaged expression across all replicates/patients.")

# OPTION B: Select a specific replicate/patient (Uncomment to use instead of A)
# target_rep = "at"  # Replace with your specific patient/replicate ID
# df_processed = df_synth.xs(target_rep, level="patient")
# print(f"  -> Selected data exclusively for: {target_rep}")
# ---------------------------------------------------------

# Extract just the target cell type (assuming cell type is now the 'celltype' level)
df_target_ct = df_processed.xs(target_ct, level="celltype")

# Setup BORF Inputs
labels = np.array([target_ct] * len(df_target_ct))
genes_vec = df_target_ct.index.get_level_values("gene").values
X_raw = df_target_ct.values[:, np.newaxis, :]

print(f"\n2. Running BORF Transformation for {target_ct}...")

# Initialize Builder (Keeping it simple, letting fast_borf handle the pipeline)
builder = BorfBuilder(
    n_jobs=-1, 
    configs=CUSTOM_CONFIG_A3_NO_DILATION_WINDOW_SIZE_3,
    pipeline_objects=[
        (ReshapeTo2D, dict(keep_unraveled_index=True)), 
        (ZeroColumnsRemover, dict(axis=0)), 
        (ToScipySparse, dict())
    ]
)

# Fit and transform
borf_model = builder.build(X_raw)
X_dense = borf_model.fit_transform(X_raw).toarray()

# X_dense to 1, if different from 0, else 0 (binary presence/absence of pattern)
X_dense = (X_dense != 0).astype(int)

print(f"BORF Matrix Shape: {X_dense.shape} (Genes x Patterns)")

# 1. Prepare your inputs
is_target_ct = (labels == target_ct)
X_dense_ct = X_dense[is_target_ct]
genes_ct = genes_vec[is_target_ct]
df_X = pd.get_dummies(genes_ct, dtype=int)

# 1. Define the generalized worker function
def process_pattern_model(idx, X_dense_ct_slice, df_X, model_name, model_template):
    y_counts = X_dense_ct_slice
    start_time = time.time()
    
    # Clone and fit the model
    model = clone(model_template)
    model.fit(df_X, y_counts)
    
    # Generate predictions for metrics
    preds = model.predict(df_X)
    
    acc_score = accuracy_score(y_counts, preds)
    prec_score = precision_score(y_counts, preds, zero_division=0)
    rec_score = recall_score(y_counts, preds, zero_division=0)
    f1_val = f1_score(y_counts, preds, zero_division=0)
    
    # Extract Feature Importance via SHAP
    if model_name in ["Decision Tree", "Random Forest", "Extra Trees"]:
        # fasttreeshap is highly optimized for sklearn tree ensembles
        explainer = fasttreeshap.TreeExplainer(model)
        shap_values_raw = explainer.shap_values(df_X, check_additivity=False)
        
    elif model_name == "LightGBM":
        # Official SHAP handles LightGBM serialization properly
        explainer = shap.TreeExplainer(model)
        shap_values_raw = explainer.shap_values(df_X, check_additivity=False)
        
    elif model_name in ["Ridge Classifier", "Logistic Regression"]:
        # Linear models use LinearExplainer
        explainer = shap.LinearExplainer(model, df_X)
        shap_values_raw = explainer.shap_values(df_X)
        
    # Safely parse SHAP formats (sklearn vs lightgbm vs linear)
    if isinstance(shap_values_raw, list):
        # SKLearn tree models return a list: [class_0, class_1]
        shap_values_class_1 = shap_values_raw[1] 
    elif shap_values_raw.ndim == 3:
        # Catch-all for 3D arrays
        shap_values_class_1 = shap_values_raw[:, :, 1]
    else:
        # LightGBM and Linear models return a single 2D array for binary classification
        shap_values_class_1 = shap_values_raw
        
    # Calculate mean positive SHAP
    importance_array = np.maximum(shap_values_class_1, 0).mean(axis=0)
        
    time_taken = time.time() - start_time
    
    return idx, y_counts.sum(), acc_score, prec_score, rec_score, f1_val, importance_array, time_taken

# 2. Setup Models (n_jobs=1 to avoid joblib parallel clashes)
models = {
    "Decision Tree": DecisionTreeClassifier(max_depth=20, random_state=42),
    "Random Forest": RandomForestClassifier(n_estimators=50, random_state=42, max_depth=None, n_jobs=1),
    "Extra Trees": ExtraTreesClassifier(n_estimators=50, random_state=42, max_depth=None, n_jobs=1),
    "LightGBM": LGBMClassifier(n_estimators=50, random_state=42, n_jobs=1, verbose=-1),
    "Ridge Classifier": RidgeClassifier(random_state=42),
    "Logistic Regression": LogisticRegression(penalty='l1', solver='liblinear', random_state=42)
}

TOP_N_GENES = 30
aggregated_results = []
gene_importances_dict = {}

print("Running parallel pattern evaluation across models...\n")

# 3. Iterate over models and run parallel processing
for model_name, model_template in models.items():
    print(f"--- Evaluating: {model_name} ---")
    
    # Execute in parallel across patterns
    results = Parallel(n_jobs=-1)(
        delayed(process_pattern_model)(idx, X_dense_ct[:, idx], df_X, model_name, model_template) 
        for idx in range(X_dense_ct.shape[1])
    )
    
    # Unpack results for the current model
    pattern_accuracies = []
    pattern_precisions = []
    pattern_recalls = []
    pattern_f1s = []
    total_time = 0
    importance_dict = {}
    
    for idx, total_count, acc, prec, rec, f1, importance_array, time_taken in results:
        pattern_accuracies.append(acc)
        pattern_precisions.append(prec)
        pattern_recalls.append(rec)
        pattern_f1s.append(f1)
        total_time += time_taken
        importance_dict[idx] = importance_array
        
    # Create DataFrame of importances (Genes x Patterns)
    importance_df = pd.DataFrame(importance_dict, index=df_X.columns)
    
    # Overall Importance: Sum across all patterns
    overall_gene_importance = importance_df.sum(axis=1).values
    gene_importances_dict[model_name] = overall_gene_importance
    
    # Extract Top Genes
    top_indices = np.argsort(overall_gene_importance)[::-1][:TOP_N_GENES]
    top_genes = df_X.columns[top_indices].tolist() 
    
    avg_acc = np.mean(pattern_accuracies)
    avg_prec = np.mean(pattern_precisions)
    avg_rec = np.mean(pattern_recalls)
    avg_f1 = np.mean(pattern_f1s)
    
    aggregated_results.append({
        "Model": model_name,
        "Mean Acc": f"{avg_acc:.4f}",
        "Mean Prec": f"{avg_prec:.4f}",
        "Mean Rec": f"{avg_rec:.4f}",
        "Mean F1": f"{avg_f1:.4f}",
        "Total Worker Time (s)": f"{total_time:.4f}",
        f"Top {TOP_N_GENES} Genes": ", ".join(top_genes)
    })
    print(f"Completed {model_name} | F1: {avg_f1:.4f} | Acc: {avg_acc:.4f}\n")

# 4. Display Ablation Results
results_df = pd.DataFrame(aggregated_results)
print("=== Ablation Study Results ===")
display(results_df)

# 5. Kendall-Tau Rank Comparison & Intersections
model_names = list(models.keys())
kendall_df = pd.DataFrame(index=model_names, columns=model_names, dtype=float)
intersection_df = pd.DataFrame(index=model_names, columns=model_names, dtype=int)

# Extract Top N Sets for easy intersection mapping
top_n_sets = {}
for m in model_names:
    top_idx = np.argsort(gene_importances_dict[m])[::-1][:TOP_N_GENES]
    top_n_sets[m] = set(df_X.columns[top_idx])

for m1, m2 in itertools.combinations(model_names, 2):
    # Kendall Tau
    tau, p_val = kendalltau(gene_importances_dict[m1], gene_importances_dict[m2])
    kendall_df.loc[m1, m2] = tau
    kendall_df.loc[m2, m1] = tau
    
    # Intersection
    shared_genes = top_n_sets[m1].intersection(top_n_sets[m2])
    intersection_df.loc[m1, m2] = len(shared_genes)
    intersection_df.loc[m2, m1] = len(shared_genes)

# Fill diagonals
np.fill_diagonal(kendall_df.values, 1.0)
np.fill_diagonal(intersection_df.values, TOP_N_GENES)

kendall_df.to_csv(os.path.join(output_folder, "kendall_tau_matrix.csv"))
intersection_df.to_csv(os.path.join(output_folder, "top_gene_intersection_matrix.csv"))