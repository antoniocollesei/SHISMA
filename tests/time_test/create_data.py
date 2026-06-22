import numpy as np
import pandas as pd
import os
import networkx as nx
import random
import csv

def generate_shisma_data(n_genes=200, n_planted=0, n_timepoints=6, n_pats=5, n_celltypes=10):
    """
    Generates completely random, independent expression data and unlinked networks
    to act as a strict null background dataset.
    """
    np.random.seed(42)
    
    cell_types = np.array([f"CT{i}" for i in range(n_celltypes)])
    target_ct = "CT0"
    padding = len(str(n_genes - 1))
    genes = np.array([f"G{i:0{padding}d}" for i in range(n_genes)])
    
    # ==========================================
    # 1. PURELY RANDOM PPI NETWORK GENERATION
    # ==========================================
    # Generates a random graph over all genes without any pre-planted structural cliques
    edges = []
    for _ in range(n_genes * 2): 
        g1, g2 = np.random.choice(genes, 2)
        if g1 != g2: 
            edges.append([g1, g2])
    ppi_df = pd.DataFrame(edges, columns=["gene1", "gene2"]).drop_duplicates()

    # ==========================================
    # 2. COMPLETELY RANDOM & INDEPENDENT EXPRESSION
    # ==========================================
    print("Generating completely random, independent expression values...")
    n_samples = n_genes * n_celltypes * n_pats
    
    # Create aligned metadata grids to generate rows
    mesh_genes, mesh_cts, mesh_pats = np.meshgrid(np.arange(n_genes), np.arange(n_celltypes), np.arange(n_pats), indexing='ij')
    flat_gene_idx = mesh_genes.flatten()
    flat_ct_idx = mesh_cts.flatten()
    flat_pats = mesh_pats.flatten()
    
    # Generate entirely independent values uniformly distributed across simulated log1p space (e.g., 0 to 4).
    # This destroys all temporal trends, covariance, and celltype specificity.
    normalized_expression = np.random.uniform(0.0, 4.0, size=(n_samples, n_timepoints))
    
    # ==========================================
    # 3. DATAFRAME FORMATTING
    # ==========================================
    df_synth = pd.DataFrame(normalized_expression, columns=[f"T{t}" for t in range(n_timepoints)])
    
    flat_genes = genes[flat_gene_idx]
    flat_cts = cell_types[flat_ct_idx]
    
    df_synth.index = [f"{g}_pat{r}_{ct}" for g, ct, r in zip(flat_genes, flat_cts, flat_pats)]
    
    return df_synth, ppi_df, [], target_ct

if __name__ == "__main__":
    sample_sizes = [100, 500, 1000, 5000]
    timepoints_list = [6, 10, 15]
    patients_list = [5, 10, 15]

    os.makedirs("data", exist_ok=True)

    for n_samples in sample_sizes:
        for n_timepoints in timepoints_list:
            for n_pats in patients_list: 
                print(f"Generating: genes={n_samples}, timepoints={n_timepoints}, patients={n_pats}")
                
                prefix = f"samples{n_samples}_time{n_timepoints}_pats{n_pats}"
                n_planted = n_samples / 10

                # Forced n_planted=0 here to match the strict null constraint
                df_synth, ppi_synth, target_genes, target_ct = generate_shisma_data(
                    n_genes=n_samples, 
                    n_planted=n_planted, 
                    n_timepoints=n_timepoints, 
                    n_pats=n_pats,        
                    n_celltypes=5
                )
                
                df_synth.to_csv(f"data/{prefix}_expression.csv")
                ppi_synth.to_csv(f"data/{prefix}_ppi.csv", index=False)