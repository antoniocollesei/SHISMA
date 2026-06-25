import numpy as np
import pandas as pd

def load_real_dataset(filepath, target_ct):
    """Loads and standardizes expression matrix indexing/columns exactly per notebook."""
    print(f"Loading real dataset from {filepath}...")
    df_real = pd.read_csv(filepath, index_col=0)
    df_real.columns = [f"T{i}" for i in range(len(df_real.columns))]
    df_real = df_real.fillna(0.0)
    print(f"  -> Loaded {df_real.shape[0]} samples across {df_real.shape[1]} timepoints.")
    print(f"  -> Target cell type set to: {target_ct}")
    return df_real, target_ct

def generate_shisma_data(n_genes=300, n_planted=30, n_timepoints=6, n_pats=6, n_celltypes=5, seed=42):
    """Generates synthetic Negative Binomial data with exact planted parameters from notebook."""
    np.random.seed(seed)
    cell_types = np.array([f"CT{i}" for i in range(n_celltypes)])
    target_ct = "CT0"
    padding = len(str(n_genes - 1))
    genes = np.array([f"G{i:0{padding}d}" for i in range(n_genes)])
    planted_genes = genes[:n_planted]
    
    # PPI network generation
    edges = [[planted_genes[i], planted_genes[j]] for i in range(n_planted) for j in range(i+1, n_planted)]
    for _ in range(n_genes * 2): 
        g1, g2 = np.random.choice(genes, 2)
        if g1 != g2: 
            edges.append([g1, g2])
    ppi_df = pd.DataFrame(edges, columns=["gene1", "gene2"]).drop_duplicates()

    # Background Expression Profile
    print("Generating realistic Negative Binomial gene expression...")
    n_samples = n_genes * n_celltypes * n_pats
    mesh_genes, mesh_cts, mesh_pats = np.meshgrid(np.arange(n_genes), np.arange(n_celltypes), np.arange(n_pats), indexing='ij')
    flat_gene_idx = mesh_genes.flatten()
    flat_ct_idx = mesh_cts.flatten()
    flat_pats = mesh_pats.flatten()
    
    gene_baselines = np.random.gamma(shape=2.0, scale=1.5, size=n_genes)
    ct_multipliers = np.random.lognormal(mean=0.0, sigma=0.5, size=(n_genes, n_celltypes))
    base_means = gene_baselines[flat_gene_idx] * ct_multipliers[flat_gene_idx, flat_ct_idx]
    
    mu_matrix = np.tile(base_means[:, np.newaxis], (1, n_timepoints))
    biological_noise = np.random.lognormal(mean=0.0, sigma=0.1, size=(n_samples, n_timepoints))
    mu_matrix = mu_matrix * biological_noise
    
    raw_counts = np.random.poisson(lam=mu_matrix)
    dropout_probs = np.exp(-0.1 * mu_matrix)
    dropouts = np.random.binomial(n=1, p=dropout_probs)
    raw_counts = raw_counts * (1 - dropouts)
    normalized_expression = np.log1p(raw_counts)
    
    # Plant target cell-type signal injection
    is_target_ct = (cell_types[flat_ct_idx] == target_ct)
    is_planted_gene = np.isin(genes[flat_gene_idx], planted_genes)
    plant_mask = is_target_ct & is_planted_gene
    normalized_expression[plant_mask, 3] += 2.0
    normalized_expression[plant_mask, 4] -= 1.5 
    normalized_expression[plant_mask, 5] += 2.0
    
    df_synth = pd.DataFrame(normalized_expression, columns=[f"T{t}" for t in range(n_timepoints)])
    df_synth.index = [f"{g}_pat{r}_{ct}" for g, ct, r in zip(genes[flat_gene_idx], cell_types[flat_ct_idx], flat_pats)]
    
    return df_synth, ppi_df, list(planted_genes), target_ct