import pandas as pd
import numpy as np
import networkx as nx
import warnings
from joblib import Parallel, delayed
from sklearn.tree import DecisionTreeClassifier
import fasttreeshap
from statsmodels.stats.multitest import multipletests
from scipy.stats import mannwhitneyu

# Import necessary modules from your fast_borf workspace dependencies
from fast_borf import BorfBuilder
from fast_borf.pipeline.zero_columns_remover import ZeroColumnsRemover
from fast_borf.pipeline.reshaper import ReshapeTo2D
from fast_borf.pipeline.to_scipy import ToScipySparse
from .parameters import BORF_CONFIG

def _process_single_pattern(idx, X_dense_ct_slice, df_X):
    """Worker function optimized to safely handle zero-variance noise and variable SHAP shapes."""
    y_counts = X_dense_ct_slice
    
    # Catch internal numpy/scikit-learn scalar division warnings from pure random noise split points
    with warnings.catch_warnings(), np.errstate(divide='ignore', invalid='ignore'):
        warnings.simplefilter("ignore")
        
        model = DecisionTreeClassifier(max_depth=20, random_state=42)
        model.fit(df_X, y_counts)
        
        explainer = fasttreeshap.TreeExplainer(model)
        shap_values_raw = explainer.shap_values(df_X, check_additivity=False)
        
        # DEFENSIVE DYNAMIC SHAPE CHECK: Safely extract class 1 values regardless of dimensions
        if isinstance(shap_values_raw, list):
            shap_values_class_1 = shap_values_raw[1] if len(shap_values_raw) > 1 else shap_values_raw[0]
        elif hasattr(shap_values_raw, "ndim") and shap_values_raw.ndim == 3:
            shap_values_class_1 = shap_values_raw[:, :, 1]
        else:
            # It is already a 2D array representing the target attribution space
            shap_values_class_1 = shap_values_raw
        
        mean_maxpos_shap = np.maximum(shap_values_class_1, 0).mean(axis=0)
        
    return idx, mean_maxpos_shap

def compute_borf_and_shap(df_synth, ppi_synth, target_ct, n_jobs=-1):
    """Transforms raw/normalized expression data using the exact notebook preprocessing and BoRF pipelines."""
    print("1. Filtering for PPI genes & Handling Replicates...")
    ppi_genes = np.unique(ppi_synth[['gene1', 'gene2']].values.astype(str))

    if not isinstance(df_synth.index, pd.MultiIndex):
        index_dims = df_synth.index.str.rsplit('_', n=2, expand=True)
        index_dims.names = ["gene", "patient", "celltype"] 
        df_synth.index = index_dims

    df_synth = df_synth[df_synth.index.get_level_values("gene").isin(ppi_genes)]
    print(f"  -> Reduced to {len(df_synth.index.get_level_values('gene').unique())} valid PPI genes.")

    df_target_ct = df_synth.xs(target_ct, level="celltype")

    labels = np.array([target_ct] * len(df_target_ct))
    genes_vec = df_target_ct.index.get_level_values("gene").values
    X_raw = df_target_ct.values[:, np.newaxis, :]

    print(f"\n2. Running BORF Transformation for {target_ct}...")
    builder = BorfBuilder(
        n_jobs=n_jobs, 
        configs=BORF_CONFIG,
        pipeline_objects=[
            (ReshapeTo2D, dict(keep_unraveled_index=True)), 
            (ZeroColumnsRemover, dict(axis=0)), 
            (ToScipySparse, dict())
        ]
    )

    borf_model = builder.build(X_raw)
    X_dense = borf_model.fit_transform(X_raw).toarray()
    X_dense = (X_dense != 0).astype(int)
    print(f"BORF Matrix Shape: {X_dense.shape} (Genes x Patterns)")

    is_target_ct = (labels == target_ct)
    X_dense_ct = X_dense[is_target_ct]
    genes_ct = genes_vec[is_target_ct]
    df_X = pd.get_dummies(genes_ct, dtype=int)

    print("Running Decision Trees and SHAP extraction in parallel via joblib...")
    results = Parallel(n_jobs=n_jobs)(
        delayed(_process_single_pattern)(idx, X_dense_ct[:, idx], df_X) 
        for idx in range(X_dense_ct.shape[1])
    )

    shap_dict = {idx: array for idx, array in results}
    shap_df = pd.DataFrame(shap_dict, index=df_X.columns)
    print("SHAP feature attribution complete!")
    return shap_df, borf_model, X_raw

def resolve_overlapping_subnetworks(df, overlap_thresh=0.5, strategy='largest', mht='fdr'):
    """Removes redundant subnetwork modules based on the Overlap Coefficient (Cell 13)."""
    if df.empty: return df
    resolved_rows = []
    
    for pattern, group in df.groupby('Pattern'):
        if strategy == 'largest':
            group = group.sort_values(by=['Size', 'Q_Value_' + mht.upper()], ascending=[False, True])
        elif strategy == 'significant':
            group = group.sort_values(by=['Q_Value_' + mht.upper(), 'Size'], ascending=[True, False])
        elif strategy == 'optimal_core':
            group = group.sort_values(by=['Q_Value_' + mht.upper(), 'Size'], ascending=[True, True])
            
        accepted_clusters = []
        for _, row in group.iterrows():
            current_genes = set(row['Genes'].split(', '))
            is_redundant = False
            
            for accepted_row in accepted_clusters:
                accepted_genes = set(accepted_row['Genes'].split(', '))
                intersection = len(current_genes.intersection(accepted_genes))
                min_size = min(len(current_genes), len(accepted_genes))
                
                if (intersection / min_size) >= overlap_thresh:
                    is_redundant = True
                    break
            
            if not is_redundant:
                accepted_clusters.append(row)
        resolved_rows.extend(accepted_clusters)
        
    return pd.DataFrame(resolved_rows)

def filter_by_cell_specificity(resolved_df, cell_expr_df, cell_metadata, target_cell_type, p_val_thresh=0.05, min_fc=1.2):
    """Filters discovered subnetworks using a Mann-Whitney U validation check (Cell 13)."""
    if resolved_df.empty: return resolved_df
    specific_subnetworks = []
    
    aligned_cells = cell_expr_df.index.intersection(cell_metadata.index)
    cell_expr_df = cell_expr_df.loc[aligned_cells]
    cell_metadata = cell_metadata.loc[aligned_cells]
    
    target_mask = (cell_metadata == target_cell_type)
    bg_mask = ~target_mask

    for _, row in resolved_df.iterrows():
        subnetwork_genes = row['Genes'].split(', ')
        valid_genes = [g for g in subnetwork_genes if g in cell_expr_df.columns]
        
        if not valid_genes: continue
            
        cell_activity_scores = cell_expr_df[valid_genes].mean(axis=1)
        target_scores = cell_activity_scores[target_mask]
        bg_scores = cell_activity_scores[bg_mask]
        
        target_mean = target_scores.mean()
        bg_mean = bg_scores.mean()
        fold_change = (target_mean + 1e-9) / (bg_mean + 1e-9)
        
        if len(target_scores) < 3 or len(bg_scores) < 3:
            p_val = 0.0
        else:
            try:
                stat, p_val = mannwhitneyu(target_scores, bg_scores, alternative='greater')
            except ValueError:
                p_val = 1.0 
            
        if p_val <= p_val_thresh and fold_change >= min_fc:
            row_dict = row.to_dict()
            row_dict['Target_Mean'] = round(target_mean, 4)
            row_dict['Background_Mean'] = round(bg_mean, 4)
            row_dict['Fold_Change'] = round(fold_change, 3)
            row_dict['Specificity_P_Val'] = p_val
            specific_subnetworks.append(row_dict)
            
    final_specific_df = pd.DataFrame(specific_subnetworks)
    if not final_specific_df.empty:
        final_specific_df = final_specific_df.sort_values(by='Fold_Change', ascending=False)
    return final_specific_df

def run_shisma_pipeline(df_synth, ppi_synth, target_ct, beta=0.1, min_size=3, max_size=30, 
                        n_perms=1000, thresholds=(50, 70, 80, 90, 95, 99), j_eps=0.01, alpha=0.05, mht='fdr', n_jobs=-1, plot=False, plot_dir=None):
    """Runs the complete, synchronized end-to-end SHISMA processing configuration."""
    
    if not isinstance(df_synth.index, pd.MultiIndex):
        index_dims = df_synth.index.str.rsplit('_', n=2, expand=True)
        index_dims.names = ["gene", "patient", "celltype"] 
        df_synth.index = index_dims

    df_backup = df_synth.copy()
    
    shap_df, borf_model, X_raw = compute_borf_and_shap(df_synth, ppi_synth, target_ct, n_jobs=n_jobs)

    G = nx.from_pandas_edgelist(ppi_synth, 'gene1', 'gene2')
    nodes = list(G.nodes())
    n_nodes = len(nodes)

    jaccard_preds = nx.jaccard_coefficient(G, G.edges())
    for u, v, j_score in jaccard_preds:
        G[u][v]['density_weight'] = j_score + j_eps

    A = nx.to_numpy_array(G, nodelist=nodes, weight='density_weight')
    col_sums = A.sum(axis=0)
    col_sums[col_sums == 0] = 1 
    W = A / col_sums

    I = np.eye(n_nodes)
    P = beta * np.linalg.inv(I - (1 - beta) * W)
    A_binary = nx.to_numpy_array(G, nodelist=nodes, weight=None)

    significant_subnetworks = []

    for idx in range(shap_df.shape[1]):
        pattern_name = f"PATTERN_{idx}"
        gene_scores = shap_df.iloc[:, idx].copy().abs()
        gene_scores.index = gene_scores.index.astype(str)

        f = np.zeros(n_nodes)
        for i, node in enumerate(nodes):
            if str(node) in gene_scores.index:
                f[i] = gene_scores[str(node)]

        s_obs = np.dot(P, f)
        S_edges_obs = A_binary * np.outer(s_obs, s_obs)
        positive_edges = S_edges_obs[S_edges_obs > 0]
        
        if len(positive_edges) == 0: continue

        candidate_clusters = set()
        for pct in thresholds:
            delta = np.percentile(positive_edges, pct)
            S_sparse = np.where(S_edges_obs >= delta, S_edges_obs, 0)
            G_thresh = nx.from_numpy_array(S_sparse, create_using=nx.Graph)
            
            for comp in nx.connected_components(G_thresh):
                if min_size <= len(comp) <= max_size:
                    candidate_clusters.add(tuple(sorted([nodes[i] for i in comp])))

        if not candidate_clusters: continue

        f_perms = np.column_stack([np.random.permutation(f) for _ in range(n_perms)])
        s_null_matrix = np.dot(P, f_perms)

        cluster_data = []
        p_values = []
        for comp in candidate_clusters:
            comp_indices = [nodes.index(g) for g in comp]
            obs_score = np.median(s_obs[comp_indices])
            null_scores = np.median(s_null_matrix[comp_indices, :], axis=0)
            
            beats_true = np.sum(null_scores >= obs_score)
            emp_p = (beats_true + 1) / (n_perms + 1)
            
            p_values.append(emp_p)
            cluster_data.append({'Size': len(comp), 'Median_Score': obs_score, 'Genes': comp})

        if mht == 'fdr':
            reject, q_values, _, _ = multipletests(p_values, alpha=alpha, method='fdr_bh')
        elif mht == 'bonferroni':
             reject, q_values, _, _ = multipletests(p_values, alpha=alpha, method='bonferroni')
        else:
            raise ValueError(f"Unsupported multiple hypothesis testing method: {mht}")

        for i, cluster in enumerate(cluster_data):
            if reject[i]:
                significant_subnetworks.append({
                    'Pattern': pattern_name,
                    'Size': cluster['Size'],
                    'Median_Score': cluster['Median_Score'],
                    'P_Value': p_values[i],
                    'Q_Value_' + mht.upper(): q_values[i],
                    'Genes': ", ".join(cluster['Genes'])
                })

    final_results_df = pd.DataFrame(significant_subnetworks) if significant_subnetworks else pd.DataFrame(columns=['Pattern', 'Size', 'Median_Score', 'P_Value', 'Q_Value_' + mht.upper(), 'Genes'])

    print("\nApplying Post-Processing Overlap Filtering...")
    filtered_results_df = resolve_overlapping_subnetworks(final_results_df, overlap_thresh=0.5, strategy='largest', mht=mht)

    print("Formatting metadata objects for cell specificity checks...")
    mean_expr = df_backup.mean(axis=1)
    df_temp = pd.DataFrame({'expr': mean_expr})
    df_temp['Gene'] = df_temp.index.get_level_values("gene")
    df_temp['Cell'] = df_temp.index.get_level_values("celltype") + '_' + df_temp.index.get_level_values("patient")
    df_temp['CellType'] = df_temp.index.get_level_values("celltype")

    cell_expr_matrix = df_temp.pivot(index='Cell', columns='Gene', values='expr')
    cell_types_series = df_temp.drop_duplicates('Cell').set_index('Cell')['CellType']

    print("Running Mann-Whitney U cell specificity filters...")
    cell_specific_df = filter_by_cell_specificity(
        resolved_df=filtered_results_df, 
        cell_expr_df=cell_expr_matrix, 
        cell_metadata=cell_types_series, 
        target_cell_type=target_ct,
        p_val_thresh=0.05,
        min_fc=1.2
    )

    # Filter final results to include only the columns of interest
    cell_specific_df = cell_specific_df[['Pattern', 'Size', 'Q_Value_' + mht.upper(), 'Genes', 'Specificity_P_Val']]

    # Rename columns for clarity
    cell_specific_df.rename(columns={
        'Q_Value_' + mht.upper(): 'Pvalue_' + mht.upper(),
        'Specificity_P_Val': target_ct + '_Specificity_Pvalue'
    }, 
    inplace=True)
    
    if plot:
        from .plotting import plot_dynamics
        if plot_dir is None:
            import os
            plot_dir = os.path.join("plots", target_ct)
        print(f"Generating dynamics plots in {plot_dir}...")
        plot_dynamics(cell_specific_df, df_backup, borf_model, X_raw, target_ct, plot_dir)
        
    return cell_specific_df