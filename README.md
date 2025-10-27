# SHISMA
**SHISMA: SHape-driven Inference of Significant celltype-specific Subnetworks from tiMe series single-cell trAnscriptomics**

Authors: Antonio Collesei, Pierangela Palmerini, Emilia Vigolo, Francesco Spinnato 

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
We strongly encourage the installation of requirements via conda/mamba. For this we make the environment .yml file available:
```bash
conda env create -f conda/SHISMA_environment.yml
conda activate SHISMA_env
```

## Usage
The dataset with time series observations must be a tab-separated csv file, with the first column representing the index names in the form of [gene]_[celltype]. Note: we will improve this in an upcoming release to allow for a more user-friendly input file.
```bash
./wrapper.sh --data <time_series_file> --ppi <ppi_file> --ct <cell_type> --out <output_name> [--nperm <num_permutations>] [--cores <num_cores>]
```

## Inputs
Time Series File (--data): Normalized expression matrix of genes across time points. \\
PPI File (--ppi): Edge list representing a biological interaction network; it must be a tab-separated file with two columns having names "gene1" and "gene2". \\
Cell Type (--ct): Specific cell type to focus the analysis on: please be sure that the exact name is present in the index of your data file.\\
Output Name (--out): Prefix for all generated output files.

## Outputs
List of significant subnetworks for the selected cell type.
Each subnetwork includes:
- A local shape (temporal pattern) signature.
- A p-value corrected for multiple testing (FWER controlled).
