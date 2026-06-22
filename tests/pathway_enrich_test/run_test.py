import pandas as pd
import numpy as np
import io
from scipy.stats import f_oneway
from statsmodels.stats.multitest import multipletests
import gseapy as gp

def run_baseline_de_and_enrichment(expr_file, ppi_file, target_ct="Bcells", pval_thresh=0.05, ppi_score_thresh=0.700):
    print(f"Loading expression data from {expr_file}...")
    df = pd.read_csv(expr_file, index_col=0)
    
    # 1. Parse the Expression Index
    index_dims = df.index.str.rsplit('_', n=2, expand=True)
    index_dims.names = ["gene", "patient", "celltype"] 
    df.index = index_dims
    
    if target_ct not in df.index.get_level_values("celltype").unique():
        print(f"Error: Cell type '{target_ct}' not found.")
        return
        
    df_ct = df.xs(target_ct, level="celltype")
    timepoints = df_ct.columns.tolist()
    
    # 2. Load and Process the STRING PPI Network
    print(f"Loading PPI network from {ppi_file}...")
    # Using delim_whitespace=True or sep='\t' depending on how your text file is saved
    ppi_df = pd.read_csv(ppi_file, delim_whitespace=True) 
    
    # Filter for high confidence interactions to match SHISMA's environment
    ppi_df = ppi_df[ppi_df['combined_score'] >= ppi_score_thresh]
    
    # Extract the valid universe of genes that exist in the network
    ppi_genes = set(ppi_df['gene1']).union(set(ppi_df['gene2']))
    print(f"PPI filtering complete: {len(ppi_genes)} unique genes retained (Score >= {ppi_score_thresh}).")
    
    # 3. Restrict Expression Data to the PPI Universe
    available_genes = set(df_ct.index.get_level_values("gene").unique())
    valid_genes = list(available_genes.intersection(ppi_genes))
    print(f"Isolated {target_ct}: {len(valid_genes)} genes are present in both the dataset and the filtered PPI.")
    
    # 4. Longitudinal Differential Expression (ANOVA)
    de_results = []
    print("Running longitudinal Differential Expression (ANOVA)...")
    for gene in valid_genes:
        gene_data = df_ct.xs(gene, level="gene")
        
        # Require at least 2 replicates to run ANOVA
        if len(gene_data) < 2:
            continue
            
        time_groups = [gene_data[tp].values for tp in timepoints]
        
        if np.all(gene_data.values == 0):
            continue
            
        stat, pval = f_oneway(*time_groups)
        
        if not np.isnan(pval):
            de_results.append({'gene': gene, 'f_stat': stat, 'p_val': pval})
            
    de_df = pd.DataFrame(de_results)
    
    if de_df.empty:
        print("No genes had sufficient variance/replicates to run DE.")
        return
        
    # 5. False Discovery Rate (FDR) Correction
    _, q_vals, _, _ = multipletests(de_df['p_val'], alpha=pval_thresh, method='fdr_bh')
    de_df['q_val_fdr'] = q_vals
    
    sig_genes_df = de_df[de_df['q_val_fdr'] < pval_thresh].sort_values('q_val_fdr')
    sig_gene_list = sig_genes_df['gene'].tolist()
    
    print(f"\nDE Analysis Complete: Found {len(sig_gene_list)} significantly time-varying genes (FDR < {pval_thresh}).")
    
    if len(sig_gene_list) == 0:
        print("No significant genes found to run enrichment on.")
        return
        
    # 6. Pathway Enrichment (Reactome 2022 via GSEApy)
    print("\nRunning Pathway Enrichment (Reactome 2022)...")
    try:
        enr = gp.enrichr(
            gene_list=sig_gene_list,
            background=list(ppi_genes), # Set background strictly to the PPI universe
            gene_sets=['Reactome_2022'],
            organism='human',
            outdir=None 
        )
        
        res_df = enr.results
        sig_pathways = res_df[res_df['Adjusted P-value'] < 0.05].sort_values('Adjusted P-value')
        
        print("\n" + "="*80)
        print(f"TOP ENRICHED PATHWAYS FOR {target_ct.upper()} (FAIR BASELINE)")
        print("="*80)
        
        if not sig_pathways.empty:
            display_cols = ['Term', 'Adjusted P-value', 'Genes']
            print(sig_pathways[display_cols].to_string(index=False))
            sig_pathways.to_csv(f"{target_ct}_longitudinal_anova_enrichment.csv", index=False)
        else:
            print("No significant pathways found (Adjusted P-value < 0.05).")
            
    except Exception as e:
        print(f"Enrichment failed: {e}")

def enrich_shisma_subnetworks(shisma_csv, ppi_file, ppi_score_thresh=0.700):
    print(f"Loading SHISMA results from {shisma_csv}...")
    shisma_df = pd.read_csv(shisma_csv)
    
    # 1. Load the exact same background universe (PPI network)
    print(f"Loading background PPI from {ppi_file}...")
    ppi_df = pd.read_csv(ppi_file, delim_whitespace=True)
    ppi_df = ppi_df[ppi_df['combined_score'] >= ppi_score_thresh]
    ppi_genes = list(set(ppi_df['gene1']).union(set(ppi_df['gene2'])))
    
    print(f"Loaded {len(shisma_df)} subnetworks. Running Reactome enrichment per module...\n")
    print("="*90)
    print(f"{'PATTERN':<12} | {'SIZE':<5} | {'TOP ENRICHED REACTOME PATHWAY (Adjusted P-val)':<60}")
    print("="*90)
    
    all_enrichment_results = []
    
    for index, row in shisma_df.iterrows():
        pattern = row['Pattern']
        size = row['Size']
        # Extract the gene list, clean up spaces
        genes = [g.strip() for g in row['Genes'].split(',')]
        
        try:
            # Run GSEApy Enrichr for this specific subnetwork
            enr = gp.enrichr(
                gene_list=genes,
                background=ppi_genes, # STRICT FAIR BASELINE
                gene_sets=['Reactome_2022'],
                organism='human',
                outdir=None
            )
            
            res_df = enr.results
            sig_pathways = res_df[res_df['Adjusted P-value'] < 0.05].sort_values('Adjusted P-value')
            
            if not sig_pathways.empty:
                # Grab the top pathway for the console summary
                top_term = sig_pathways.iloc[0]['Term']
                top_pval = sig_pathways.iloc[0]['Adjusted P-value']
                print(f"{pattern:<12} | {size:<5} | {top_term[:55]:<55} ({top_pval:.2e})")
                
                # Save all significant pathways for this subnetwork for deeper analysis
                sig_pathways = sig_pathways.copy()
                sig_pathways['Pattern'] = pattern
                sig_pathways['Subnetwork_Size'] = size
                sig_pathways['Subnetwork_Genes'] = row['Genes']
                all_enrichment_results.append(sig_pathways)
            else:
                print(f"{pattern:<12} | {size:<5} | {'No significant pathways found (<0.05)':<55} (-)")
                
        except Exception as e:
            print(f"{pattern:<12} | {size:<5} | Error during enrichment: {e}")
            
    if all_enrichment_results:
        final_enrichment_df = pd.concat(all_enrichment_results, ignore_index=True)
        final_enrichment_df.to_csv("shisma_subnetworks_enriched.csv", index=False)
        print("\nSaved full subnetwork enrichment details to 'shisma_subnetworks_enriched.csv'.")


if __name__ == "__main__":
    run_baseline_de_and_enrichment(
        expr_file="../../data/temporal_data_with_patient_ready_normalized_full_genes.csv", 
        ppi_file="../../data/string_is_0.7_ev_reactome.tsv",
        target_ct="Bcells",
        ppi_score_thresh=0.700
    )

    enrich_shisma_subnetworks(
        shisma_csv="../../results/Bcells_significant_subnetworks_fdr.csv", 
        ppi_file="../../data/string_is_0.7_ev_reactome.tsv",
        ppi_score_thresh=0.700
    )