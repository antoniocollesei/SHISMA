import os
import glob
import pandas as pd
import time

# Import the core pipeline directly from your new package
from shisma.core import run_shisma_pipeline

def main():
    data_dir = "data"
    results_dir = "results"
    
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

    # List and path to track and save execution times
    timing_records = []
    timing_csv_path = os.path.join(results_dir, "dataset_execution_times.csv")
    perm_list = [10, 100, 500, 1000]

    total_time_start = time.time()

    for expr_file in expression_files:
        # Extract the prefix (e.g., samples100_time6_pats1)
        prefix = os.path.basename(expr_file).replace("_expression.csv", "")
        ppi_file = os.path.join(data_dir, f"{prefix}_ppi.csv")

        print("="*60)
        print(f"Processing Dataset: {prefix}")
        print("="*60)

        # 1. LOAD DATA ONCE PER DATASET (Outside the permutations loop!)
        print("Loading data into memory...")
        try:
            df_expr = pd.read_csv(expr_file, index_col=0)
            ppi_df = pd.read_csv(ppi_file)
        except Exception as e:
            print(f"  ❌ Error loading data for {prefix}: {e}")
            continue

        for n_perms in perm_list:
            print(f"  -> Running {n_perms} permutations...")

            # Start the timer exactly when the algorithm begins
            start_time = time.time()

            # 2. Run the End-to-End SHISMA Pipeline
            try:
                results_df = run_shisma_pipeline(
                    df_synth=df_expr,        
                    ppi_synth=ppi_df,        
                    target_ct=target_cell_type,
                    beta=0.1,
                    min_size=3,
                    max_size=30,
                    n_perms=n_perms,  # Pass the dynamic loop variable here
                    thresholds=(50, 70, 80, 90, 95, 99),
                    alpha=0.05
                )
            except Exception as e:
                print(f"  ❌ Algorithm failed on {prefix} (Perms: {n_perms}): {e}")
                continue
            
            # Stop timer
            elapsed = time.time() - start_time
            print(f"     ✅ Finished in {elapsed:.2f} seconds.")
            
            # Log the timing for this dataset + permutation combo
            timing_records.append({
                "Dataset": prefix,
                "Num_Permutations": n_perms,
                "Elapsed_Time_Seconds": elapsed
            })
            
            # Incrementally save the timing DataFrame to disk
            pd.DataFrame(timing_records).to_csv(timing_csv_path, index=False)
        
    total_elapsed = time.time() - total_time_start
    print("="*60)
    print(f"Batch processing complete! Processed {len(expression_files)} datasets across {len(perm_list)} permutation settings in {total_elapsed:.2f} seconds.")
    print(f"Timing logs saved to: {timing_csv_path}")
    print("="*60)

if __name__ == "__main__":
    main()