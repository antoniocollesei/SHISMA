import numpy as np
import pandas as pd
import networkx as nx
import random
import csv

np.random.seed(0)

import os
import numpy as np
import pandas as pd
import networkx as nx
import random
import csv


def generate_expression_data(
    n_samples=1000,
    n_genes=100,
    n_timepoints=6,
    cell_types=("Bcell", "Tcell", "NKcell", "Monocyte", "Dendritic"),
    filename="random_data.csv",
):
    """Generate a random gene expression dataset with given parameters."""
    # Expression values
    data = np.random.rand(n_samples, n_timepoints)
    columns = [f"t{i}" for i in range(n_timepoints)]

    # Gene and cell type labeling
    genes = [f"gene{i}" for i in range(n_genes)]
    assigned_genes = [genes[i % n_genes] for i in range(n_samples)]
    assigned_cells = [cell_types[i % len(cell_types)] for i in range(n_samples)]
    index = [f"{g}_{c}" for g, c in zip(assigned_genes, assigned_cells)]

    df = pd.DataFrame(data, columns=columns, index=index)
    df.to_csv(filename)
    return df


def generate_ppi(
    n_genes=100,
    min_extra=None,
    max_extra=None,
    filename="random_ppi.csv",
):
    """Generate a connected random PPI network and save as CSV."""
    genes = [f"gene{i}" for i in range(n_genes)]
    G = nx.Graph()
    G.add_nodes_from(genes)

    # Ensure connectivity with spanning tree
    nodes = genes[:]
    random.shuffle(nodes)
    for i in range(len(nodes) - 1):
        G.add_edge(nodes[i], nodes[i + 1])

    # Add more edges
    if min_extra is None:
        min_extra = n_genes
    if max_extra is None:
        max_extra = 2 * n_genes

    extra_edges = random.randint(min_extra, max_extra)
    for _ in range(extra_edges):
        u, v = random.sample(genes, 2)
        G.add_edge(u, v)

    # Save to file
    with open(filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow(["gene1", "gene2"])
        for edge in G.edges():
            writer.writerow([edge[0], edge[1]])
    return G


def generate_dataset_ppi(
    n_samples=1000,
    n_genes=100,
    n_timepoints=6,
    prefix="dataset1",
    outdir="datasets",
):
    """Generate a dataset + PPI pair with given parameters, saved in a folder."""
    os.makedirs(outdir, exist_ok=True)

    expr_file = os.path.join(outdir, f"{prefix}.csv")
    ppi_file = os.path.join(outdir, f"{prefix}_ppi.tsv")

    df = generate_expression_data(
        n_samples=n_samples,
        n_genes=n_genes,
        n_timepoints=n_timepoints,
        filename=expr_file,
    )

    G = generate_ppi(n_genes=n_genes, filename=ppi_file)

    return df, G

if __name__ == "__main__":
    sample_sizes = [100, 500, 1000, 2000]
    timepoints_list = [4, 6, 8, 10]
    n_genes = 20

    for n_samples in sample_sizes:
        for n_timepoints in timepoints_list:
            prefix = f"samples{n_samples}_time{n_timepoints}"
            generate_dataset_ppi(
                n_samples=n_samples,
                n_genes=n_genes,
                n_timepoints=n_timepoints,
                prefix=prefix,
                outdir="datasets"
            )

