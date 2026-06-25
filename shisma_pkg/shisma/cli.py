# shisma/cli.py

import argparse
import os
import pandas as pd
import shisma

def main():
    parser = argparse.ArgumentParser(
        description="SHISMA CLI Engine: Extract Functional Networks from Time-Series Data."
    )
    parser.add_argument("--expr-matrix", type=str, help="Path to normalized gene expression CSV dataset.")
    parser.add_argument("--ppi-network", type=str, help="Path to interacting node network table file (CSV/TSV).")
    parser.add_argument("--target-ct", type=str, required=True, help="Target specific context cell-type identifier code name.")
    parser.add_argument("--output-dir", type=str, default="shisma_results", help="Directory destination to dispatch calculated files/plots.")
    parser.add_argument("--synthetic", action="store_true", help="Generate and evaluate synthetic Negative Binomial pipeline data.")
    parser.add_argument("--borf-config", type=str, default=None, help="Path to an alternative configuration file.")
    parser.add_argument("--beta", type=float, default=0.1, help="Network affinity propagation restart damping rate variable.")
    parser.add_argument("--overlap-thresh", type=float, default=0.5, help="Cluster merge filter overlap coefficient parameter.")
    parser.add_argument("--n-perms", type=int, default=1000, help="Total permutation shuffles to perform for null matrix testing.")
    parser.add_argument("--mht-correction", type=str, default="fdr", help="Multiple hypothesis testing correction method.")
    parser.add_argument("--alpha", type=float, default=0.05, help="Significance level for multiple hypothesis testing correction.")
    parser.add_argument("--n-jobs", type=int, default=-1, help="Number of parallel jobs to run for BORF and SHAP computations.")
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    if args.synthetic:
        print("Executing standard internal simulation data profile step...")
        df_expression, ppi_df, _, target_ct = shisma.generate_shisma_data()
    else:
        if not args.expr_matrix or not args.ppi_network:
            parser.error("You must explicitly supply values for parameters: --expr-matrix and --ppi-network unless using option switch: --synthetic.")
        
        df_expression, target_ct = shisma.load_real_dataset(args.expr_matrix, args.target_ct)
        sep = '\t' if args.ppi_network.endswith('.tsv') else ','
        ppi_df = pd.read_csv(args.ppi_network, sep=sep)[['gene1', 'gene2']]
        
    # Unpack expanded array objects out from transformation execution layer
    X_binary, genes_vec, borf_model, X_raw, df_target_ct = shisma.run_borf_transformation(
        df_expression, ppi_df, target_ct, config_file_path=args.borf_config, n_jobs=args.n_jobs
    )
    
    # Calculate feature importances matrix properties via parallel TreeSHAP solvers
    shap_df = shisma.compute_shap_importance(X_binary, genes_vec, target_ct, n_jobs=args.n_jobs)
    
    # Topological smoothing matrix operations kernel 
    extractor = shisma.NetworkExtractor(ppi_df, beta=args.beta)
    subnetworks_df = extractor.extract_subnetworks(shap_df, n_perms=args.n_perms, overlap_thresh=args.overlap_thresh, mht_correction=args.mht_correction, alpha=args.alpha)
    
    # Isolate non-redundant cell-type specific target feature domains 
    specific_df = shisma.filter_specificity(subnetworks_df, df_expression, target_ct)
    
    if not specific_df.empty:
        print("\n" + "="*70 + "\nFINAL DISCOVERED SPECIFIC NETWORK STRUCTURES\n" + "="*70)
        print(specific_df[['Pattern', 'Size', 'P_Value_adj', 'Target_Mean', 'Background_Mean', 'Fold_Change', 'Genes']].to_string(index=False))
        
        # Save output results table CSV file, after resetting index
        specific_df = specific_df.reset_index(drop=True)
        specific_df.insert(0, 'Subnetwork_Index', specific_df.index) # Add index column for easier tracing in plots
        specific_df.to_csv(os.path.join(args.output_dir, f"{target_ct}_significant_subnetworks.csv"), index=False)
        
        # Launch corrected side-by-side plotting pipeline engine
        print("\nLaunching visualization pipeline engine to generate side-by-side trend layouts...")
        shisma.plot_all_subnetworks(
            cell_specific_df=specific_df,
            df_target_ct=df_target_ct,
            borf_model=borf_model,
            X_raw=X_raw,
            target_ct=target_ct,
            output_dir=args.output_dir
        )
    else:
        print("\nProcess termination: No significant modules satisfied the filtering threshold profiles.")

if __name__ == "__main__":
    main()