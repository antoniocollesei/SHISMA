# SHISMA
**SHISMA: SHape-driven Inference of Significant celltype-specific Subnetworks from tiMe series single-cell trAnscriptomics**

Authors: Antonio Collesei, Francesco Spinnato, Pierangela Palmerini, Emilia Vigolo  
Conference: Proceedings of the 20th Conference on Computational Intelligence Methods for Bioinformatics & Biostatistics (CIBB 2025)

---

## Overview
**SHISMA** is a novel algorithm designed to infer statistically significant, cell type-specific subnetworks from time-series single-cell RNA-seq data.  
It combines time series pattern mining (Bag-of-Receptive-Fields) with graph algorithms on Protein-Protein Interaction (PPI) networks, ensuring high statistical rigor through Family-Wise Error Rate (FWER) correction.

## Installation
Clone the repository:
```bash
git clone https://github.com/antoniocollesei/SHISMA.git
cd SHISMA
```
We strongly encourage the installation of requirements via conda/mamba:
```bash
conda env create -f conda/environment.yml
conda activate SHISMA_env
```

## Usage

```bash
./wrapper.sh --data <time_series_file> --ppi <ppi_file> --ct <cell_type> --out <output_name> [--nperm <num_permutations>] [--cores <num_cores>]
```

## Inputs
Time Series File (--data): Normalized expression matrix of genes across time points.
PPI File (--ppi): Edge list or adjacency matrix representing a biological interaction network.
Cell Type (--ct): Specific cell type to focus the analysis on.
Output Name (--out): Prefix for all generated output files.

## Outputs
List of significant subnetworks for the selected cell type.
Each subnetwork includes:
- A local shape (temporal pattern) signature.
- A p-value corrected for multiple testing (FWER controlled).

## Working Example
```bash
./wrapper.sh --data ../tests/random_data_test/random_data.csv \
             --ppi ../tests/random_data_test/random_ppi.csv \
             --ct Bcell \
             --out results_bcells \
             --nperm 100 \
             --cores 16
```
