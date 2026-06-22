import argparse
import os
import pandas as pd
from .core import run_shisma_pipeline

def main():
    parser = argparse.ArgumentParser(description="SHISMA CLI: Fully Aligned End-to-End Execution Pipeline")
    
    parser.add_argument("--data", required=True, help="Path to normalized expression matrix")
    parser.add_argument("--ppi", required=True, help="Path to PPI network file (sep row: gene1, gene2)")
    parser.add_argument("--ct", required=True, help="Target celltype to isolate")
    parser.add_argument("--out", required=True, help="Output CSV results file")
    
    parser.add_argument("--beta", type=float, default=0.1, help="RWR restart probability factor")
    parser.add_argument("--min_size", type=int, default=3, help="Minimum constraint bounds for components")
    parser.add_argument("--max_size", type=int, default=30, help="Maximum constraint bounds for components")
    parser.add_argument("--nperm", type=int, default=1000, help="Total permutation matrices to map")
    parser.add_argument("--cores", type=int, default=-1, help="Number of CPU cores to use")
    parser.add_argument("--j_eps", type=float, default=0.01, help="Minimal Jaccard margin score baseline")
    parser.add_argument("--alpha", type=float, default=0.05, help="BH FDR Q-value cutoff constraint threshold")
    parser.add_argument("--mht", default="fdr", help="FWER correction method for multiple hypothesis testing (e.g., 'fdr', 'bonferroni')")
    parser.add_argument("--thresholds", type=int, nargs='+', default=[50, 70, 80, 90, 95, 99], help="Percentile array targets")

    args = parser.parse_args()

    print(f"Reading target matrix log file: {args.data}")
    df_expr = pd.read_csv(args.data, index_col=0)

    print(f"Reading structural connectivity list: {args.ppi}")
    ppi_df = pd.read_csv(args.ppi)
    if '\t' in open(args.ppi).readline():
        ppi_df = pd.read_csv(args.ppi, sep='\t')
    
    results_df = run_shisma_pipeline(
        df_synth=df_expr,
        ppi_synth=ppi_df,
        target_ct=args.ct,
        beta=args.beta,
        min_size=args.min_size,
        max_size=args.max_size,
        n_perms=args.nperm,
        thresholds=tuple(args.thresholds),
        j_eps=args.j_eps,
        alpha=args.alpha,
        mht=args.mht,
        n_jobs=args.cores
    )

    print("\n" + "="*70 + "\nExecution Complete! Final Significant Targets Summary:\n" + "="*70)
    if not results_df.empty:
        print(results_df.to_string(index=False))
        os.makedirs(os.path.dirname(args.out) if os.path.dirname(args.out) else '.', exist_ok=True)
        results_df.to_csv(args.out, index=True)
        print(f"\nSaved curated outputs to target destination: {args.out}")
    else:
        print("No modules passed the full pipeline filters.")

if __name__ == "__main__":
    main()