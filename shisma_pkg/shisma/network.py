import numpy as np
import pandas as pd
import networkx as nx
from statsmodels.stats.multitest import multipletests

class NetworkExtractor:
    def __init__(self, ppi_df, beta=0.1, jaccard_epsilon=0.01):
        self.beta = beta
        self.jaccard_epsilon = jaccard_epsilon
        self.G = nx.from_pandas_edgelist(ppi_df, 'gene1', 'gene2')
        self.nodes = list(self.G.nodes())
        self.n_nodes = len(self.nodes)
        
        print(f"PPI Graph: {self.n_nodes} nodes, {self.G.number_of_edges()} edges")
        self._weight_edges_by_jaccard()
        self.P = self._precompute_propagation_matrix()

    def _weight_edges_by_jaccard(self):
        print("\nWeighting edges based on local topological density (Jaccard Similarity)...")
        jaccard_preds = nx.jaccard_coefficient(self.G, self.G.edges())
        for u, v, j_score in jaccard_preds:
            self.G[u][v]['density_weight'] = j_score + self.jaccard_epsilon

    def _precompute_propagation_matrix(self):
        print("Calculating closed-form Similarity Matrix (P)...")
        A = nx.to_numpy_array(self.G, nodelist=self.nodes, weight='density_weight')
        col_sums = A.sum(axis=0)
        col_sums[col_sums == 0] = 1 
        W = A / col_sums
        I = np.eye(self.n_nodes)
        return self.beta * np.linalg.inv(I - (1 - self.beta) * W)

    def extract_subnetworks(self, shap_df, threshold_sweep=[50, 70, 80, 90, 95, 99], 
                            min_size=3, max_size=30, n_perms=1000, overlap_thresh=0.5,
                            mht_correction='fdr', alpha=0.05):
        significant_modules = []
        A_binary = nx.to_numpy_array(self.G, nodelist=self.nodes, weight=None)

        for col_idx in range(shap_df.shape[1]):
            pattern_name = f"PATTERN_{col_idx}"
            print(f"\n{'='*40}\n--- {pattern_name} ---\n{'='*40}")
            
            gene_scores = shap_df.iloc[:, col_idx].abs()
            f = np.array([gene_scores.get(str(node), 0.0) for node in self.nodes])
            
            s_obs = np.dot(self.P, f)
            S_edges_obs = A_binary * np.outer(s_obs, s_obs)
            positive_edges = S_edges_obs[S_edges_obs > 0]
            if len(positive_edges) == 0: continue

            candidate_clusters = set()
            for pct in threshold_sweep:
                delta = np.percentile(positive_edges, pct)
                G_thresh = nx.from_numpy_array(np.where(S_edges_obs >= delta, S_edges_obs, 0))
                for comp in nx.connected_components(G_thresh):
                    if min_size <= len(comp) <= max_size:
                        candidate_clusters.add(tuple(sorted([self.nodes[i] for i in comp])))

            if not candidate_clusters:
                print("  No candidate clusters formed across any threshold scale.")
                continue

            print(f"  Generated {len(candidate_clusters)} unique candidate subnetworks across scales.")
            print("  Running vectorized permutation testing...")
            
            f_perms = np.column_stack([np.random.permutation(f) for _ in range(n_perms)])
            s_null_matrix = np.dot(self.P, f_perms)

            cluster_data, p_values = [], []
            for comp in candidate_clusters:
                comp_indices = [self.nodes.index(g) for g in comp]
                obs_score = np.median(s_obs[comp_indices])
                null_scores = np.median(s_null_matrix[comp_indices, :], axis=0)
                
                emp_p = (np.sum(null_scores >= obs_score) + 1) / (n_perms + 1)
                p_values.append(emp_p)
                cluster_data.append({'Size': len(comp), 'Median_Score': obs_score, 'Genes': comp})

            if mht_correction == 'fdr':
                reject, q_values, _, _ = multipletests(p_values, alpha=alpha, method='fdr_bh')
            if mht_correction == 'bonferroni':
                reject, q_values, _, _ = multipletests(p_values, alpha=alpha, method='bonferroni')
            else:
                print(f"  Warning: Unknown multiple hypothesis correction method '{mht_correction}'. Defaulting to FDR.")
                reject, q_values, _, _ = multipletests(p_values, alpha=alpha, method='fdr_bh')

            passed_count = 0
            for i, cluster in enumerate(cluster_data):
                if reject[i]:
                    passed_count += 1
                    print(f"  ✅ Significant Subnetwork (Size {cluster['Size']}, Q-Val: {q_values[i]:.4f})")
                    significant_modules.append({
                        'Pattern': pattern_name, 'Size': cluster['Size'],
                        'Median_Score': cluster['Median_Score'], 'P_Value': p_values[i],
                        'P_Value_adj': q_values[i], 'Genes': ", ".join(cluster['Genes'])
                    })
            if passed_count == 0:
                print("  ❌ No candidate subnetworks survived FDR correction.")

        print("\n" + "="*70 + "\nFINAL SIGNIFICANT SUBNETWORKS DATAFRAME (FDR CORRECTED)\n" + "="*70)
        final_df = pd.DataFrame(significant_modules)
        return self._resolve_overlaps(final_df, overlap_thresh)

    def _resolve_overlaps(self, df, overlap_thresh):
        if df.empty: return df
        resolved_rows = []
        for pattern, group in df.groupby('Pattern'):
            group = group.sort_values(by=['P_Value_adj', 'Size'], ascending=[True, False])
            accepted_clusters = []
            for _, row in group.iterrows():
                current_genes = set(row['Genes'].split(', '))
                is_redundant = False
                for accepted_row in accepted_clusters:
                    accepted_genes = set(accepted_row['Genes'].split(', '))
                    overlap_coef = len(current_genes.intersection(accepted_genes)) / min(len(current_genes), len(accepted_genes))
                    if overlap_coef >= overlap_thresh:
                        is_redundant = True
                        break
                if not is_redundant:
                    accepted_clusters.append(row)
            resolved_rows.extend(accepted_clusters)
        return pd.DataFrame(resolved_rows)