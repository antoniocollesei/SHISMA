import pandas as pd
import numpy as np
import os

def create_subdataset(df, n_rows, n_cols, ppi, random_state=None):
    """
    Create a sub-dataset with specified rows and columns.
    Row names (index) are preserved in the output.
    """
    rng = np.random.default_rng(random_state)

    # Sample rows by index
    sampled_rows = df.sample(n=n_rows, replace=False, random_state=random_state)

    # select some genes from PPI to ensure they are included
    ppi_genes = set(ppi['gene1']).union(set(ppi['gene2']))
    # Select genes from ppi_genes
    selected_genes = rng.choice(list(ppi_genes), size=int(min(n_rows / 10, len(ppi_genes))), replace=False)

    # Define three celltypes to repeat
    celltypes = ['Bcells', 'Tcells', 'CD8cells', 'NKcells', 'Fakecells']

    if n_cols <= df.shape[1]:
        # Randomly pick a subset of columns
        chosen_cols = rng.choice(df.columns, size=n_cols, replace=False)

        # Create fake index by combining selected_genes and celltypes, repeated as needed
        fake_index = []
        for gene in selected_genes:
            for ct in celltypes:
                fake_index.append(f"{gene}_at_{ct}")
        # Repeat or trim to match n_rows
        fake_index = (fake_index * ((n_rows // len(fake_index)) + 1))[:n_rows]

        sampled_rows.index = fake_index
        return sampled_rows.loc[:, chosen_cols]
    else:
        # Oversample columns (with replacement) if n_cols > original
        chosen_cols = rng.choice(df.columns, size=n_cols, replace=True)
        expanded = pd.concat(
            [sampled_rows[col].rename(f"{col}_{i}") for i, col in enumerate(chosen_cols, 1)],
            axis=1
        )

        fake_index = []
        for gene in selected_genes:
            for ct in celltypes:
                fake_index.append(f"{gene}_at_{ct}")
        # Repeat or trim to match n_rows
        fake_index = (fake_index * ((n_rows // len(fake_index)) + 1))[:n_rows]
        expanded.index = fake_index
        return expanded
    
# Example usage
df = pd.read_csv("../../SHISMA_main/temporal_data_with_patient_ready_normalized_full_genes.csv", sep=",", index_col=0)
ppi_full = pd.read_csv("../../SHISMA_main/string_is_0.7_ev_reactome.tsv", sep="\t", header=0)

unique_genes = pd.unique(pd.concat([ppi_full['gene1'], ppi_full['gene2']]))

# sample half of the genes in ppi_full
ppi_half = ppi_full.sample(frac=0.5, random_state=42)
ppi_halfhalf = ppi_full.sample(frac=0.25, random_state=42)

# create output folder
newpath = 'datasets' 
if not os.path.exists(newpath):
    os.makedirs(newpath)

rows = [1000, 5000, 10000, 50000, 100000, 500000, 1000000]
cols = [5, 10, 15]

for r in rows:
    for c in cols:
        for ppi, ppi_name in zip([ppi_full, ppi_half, ppi_halfhalf], ['full', 'half', 'quarter']):
            sub_df = create_subdataset(df, ppi=ppi, n_rows=r, n_cols=c, random_state=42)
            sub_df.to_csv(f'{newpath}/subdataset_{ppi_name}_{r}rows_{c}cols.csv', index=True)
            ppi.to_csv(f'{newpath}/ppi_{ppi_name}_{r}rows_{c}cols.tsv', sep="\t", index=False)