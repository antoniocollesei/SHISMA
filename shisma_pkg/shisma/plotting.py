import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from fast_borf.xai.mapping import BagOfReceptiveFields

def plot_dynamics(cell_specific_df, df_target_ct, borf_model, X_raw, target_ct, output_dir):
    """
    Generates dynamic expression plots and mapper curves for significant subnetworks.
    """
    if cell_specific_df.empty:
        print("No subnetworks to plot.")
        return

    os.makedirs(output_dir, exist_ok=True)
    
    # Export mappings
    mapper = BagOfReceptiveFields(borf_model)
    mapper.build(X_raw)

    # Loop through every row (subnetwork) in the dataset
    for index, row in cell_specific_df.iterrows():
        # Extract and clean the list of genes for the current row
        top_subnetwork_genes = [g.strip() for g in row["Genes"].split(",")]
        print(f"\n--- Processing Row/Subnetwork {index} ---")
        print(f"Genes: {top_subnetwork_genes}")

        # Extract the integer ID directly from the 'Pattern' column (e.g., 'PATTERN_3' -> 3)
        match = re.search(r"\d+", str(row["Pattern"]))
        if match:
            pattern_id = int(match.group())
        else:
            print(f"Skipping: Could not parse pattern ID from Pattern column value '{row.get('Pattern', 'Missing')}'")
            continue

        # Verify the presence of genes in the MultiIndex
        gene_level_values = df_target_ct.index.get_level_values("gene")
        genes_in_data = [
            g for g in top_subnetwork_genes if g in gene_level_values
        ]

        if genes_in_data:
            # 1. Filter the original DataFrame (keeping all replicates)
            df_filtered = df_target_ct[gene_level_values.isin(genes_in_data)]

            # 2. Reset the index to convert the 'gene' level into a normal column
            df_reset = df_filtered.reset_index(level="gene")

            # 3. Transform the DataFrame from "wide" to "long" format (Tidy Data)
            time_cols = [c for c in df_reset.columns if c != "gene"]

            df_long = df_reset.melt(
                id_vars=["gene"],
                value_vars=time_cols,
                var_name="Timepoints",
                value_name="Expression",
            )

            # 4. Generate side-by-side plots
            fig, (ax1, ax2) = plt.subplots(
                1,
                2,
                figsize=(15, 6),
                gridspec_kw={"width_ratios": [3, 1]},
                constrained_layout=True,
            )

            # --- LEFT PLOT: Expression Lineplot ---
            sns.lineplot(
                data=df_long,
                x="Timepoints",
                y="Expression",
                hue="gene",
                marker="o",
                errorbar=("ci", 95), 
                ax=ax1,
            )

            # Updated title to explicitly mention the extracted Pattern ID
            ax1.set_title(
                f"Pattern {pattern_id}: Expression in {target_ct} Over Time (95% CI)"
            )
            ax1.set_xlabel("Timepoints")
            ax1.set_ylabel("Expression (Log1p)")
            
            # Dynamically calculate legend columns (up to 6 columns wide)
            col_count = min(6, len(genes_in_data))
            
            # Place legend below the plot, centered, using multiple columns
            ax1.legend(
                loc="upper center", 
                bbox_to_anchor=(0.5, -0.15), 
                title="Gene",
                ncol=col_count,
                fontsize="small",
                frameon=False 
            )
            ax1.grid(True, linestyle="--", alpha=0.5)

            # --- RIGHT PLOT: Mapper Word Array Curve ---
            try:
                # Retrieve the array from your mapper object
                word_arr = mapper[pattern_id].word_array
                x_positions = range(len(word_arr))

                # Plot the array as a curve
                ax2.plot(
                    x_positions, 
                    word_arr, 
                    marker='o', 
                    linestyle='-', 
                    color='teal', 
                    linewidth=2, 
                    markersize=6
                )
                
                ax2.set_title(f"Pattern {pattern_id} Curve")
                ax2.set_xlabel("Array Index")
                ax2.set_ylabel("Value")
                
                # Format X-axis to show discrete positions (P0, P1, P2...)
                ax2.set_xticks(x_positions)
                ax2.set_xticklabels([f"P{k}" for k in x_positions])
                
                # Prevent flat lines from hugging the top/bottom axes
                y_min, y_max = min(word_arr), max(word_arr)
                if y_min == y_max:
                    ax2.set_ylim(y_min - 1, y_max + 1)
                else:
                    padding = (y_max - y_min) * 0.1
                    ax2.set_ylim(y_min - padding, y_max + padding)
                
                ax2.grid(True, linestyle="--", alpha=0.5)

            except Exception as e:
                ax2.text(
                    0.5,
                    0.5,
                    f"Error loading mapper[{pattern_id}]:\n{str(e)}",
                    ha="center",
                    va="center",
                )
                ax2.axis("off")

            # --- SAVE ---
            # Save the figure to the output directory
            filename = os.path.join(output_dir, f"Subnetwork_{index}_Pattern_{pattern_id}.png")
            plt.savefig(filename, dpi=300, bbox_inches="tight")
            print(f"Saved plot to: {filename}")
            
            # Close the figure to avoid holding memory
            plt.close(fig)
