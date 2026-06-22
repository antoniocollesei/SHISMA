import os
import glob
import pandas as pd
import time

# Import the core pipeline directly from your new package
from shisma.core import run_shisma_pipeline

def main():
    data_dir = "random_data"
    results_dir = "random_results"
    
    # Create the results directory if it doesn't exist
    os.makedirs(results_dir, exist_ok=True)

    # Find all generated expression datasets
    expression_files = glob.glob(os.path.join(data_dir, "*_expression.csv"))
    
    if not expression_files:
        print(f"No expression datasets found in '{data_dir}'. Please run the generation script first.")
        return

    print(f"Found {len(expression_files)} synthetic datasets. Starting batch processing...\n")

    # Define the target cell type you set in your generation script
    target_cell_type = "CT0"

    total_time_start = time.time()

    for expr_file in expression_files:
        # Extract the prefix (e.g., samples100_time6_pats1)
        prefix = os.path.basename(expr_file).replace("_expression.csv", "")
        ppi_file = os.path.join(data_dir, f"{prefix}_ppi.csv")
        out_file = os.path.join(results_dir, f"{prefix}_significant_subnetworks.csv")

        print("="*60)
        print(f"Processing Dataset: {prefix}")
        print("="*60)

        # 1. Load the raw synthetic data
        try:
            df_expr = pd.read_csv(expr_file, index_col=0)
            ppi_df = pd.read_csv(ppi_file)
        except Exception as e:
            print(f"  ❌ Error loading data for {prefix}: {e}")
            continue

        start_time = time.time()

        # 2. Run the End-to-End SHISMA Pipeline
        print(f"Running pipeline for target cell type: {target_cell_type}...")
        try:
            results_df = run_shisma_pipeline(
                df_synth=df_expr,        
                ppi_synth=ppi_df,       
                target_ct=target_cell_type,
                beta=0.1,
                min_size=3,
                max_size=30,
                n_perms=1000, 
                thresholds=(50, 70, 80, 90, 95, 99),
                alpha=0.05
            )
        except Exception as e:
            print(f"  ❌ Algorithm failed on {prefix}: {e}")
            continue

        # 3. Save Results
        results_df.to_csv(out_file, index=False)
        
        elapsed = time.time() - start_time
        
        if not results_df.empty:
            print(f"\n  Found {len(results_df)} significant subnetworks in {elapsed:.2f} seconds.")
            print(f"  Saved to {out_file}\n")
        else:
            print(f"\n  No significant subnetworks found in {elapsed:.2f} seconds.")
            print(f"  Saved empty log to {out_file}\n")

    total_elapsed = time.time() - total_time_start
    print("="*60)
    print(f"Batch processing complete! Processed {len(expression_files)} datasets in {total_elapsed:.2f} seconds.")
    print("="*60)

if __name__ == "__main__":
    main()