# shisma/visualization.py

import os
import re
import matplotlib.pyplot as plt
import seaborn as sns
from fast_borf.xai.mapping import BagOfReceptiveFields

def plot_all_subnetworks(cell_specific_df, df_target_ct, borf_model, X_raw, target_ct, output_dir):
    """Generates the dual-panel line plots matching the notebook layout, preserving row tracing indices."""
    current_output_dir = os.path.join(output_dir, target_ct)
    os.makedirs(current_output_dir, exist_ok=True)

    # Instantiate the mapper object using internal model layers
    mapper = BagOfReceptiveFields(borf_model)
    mapper.build(X_raw)

    if cell_specific_df.empty:
        print("The cell_specific_df is empty. No subnetworks to plot.")
        return

    # Loop through every row preserving the exact dataframe index for tracing links
    for index, row in cell_specific_df.iterrows():
        top_subnetwork_genes = [g.strip() for g in row["Genes"].split(",")]
        print(f"\n--- Processing Row/Subnetwork {index} ---")
        print(f"Genes: {top_subnetwork_genes}")

        # Parse out pattern digits
        match = re.search(r"\d+", str(row["Pattern"]))
        if match:
            pattern_id = int(match.group())
        else:
            print(f"Skipping: Could not parse pattern ID from column: '{row.get('Pattern', 'Missing')}'")
            continue

        gene_level_values = df_target_ct.index.get_level_values("gene")
        genes_in_data = [g for g in top_subnetwork_genes if g in gene_level_values]

        if genes_in_data:
            # 1. Slice target DataFrame preserving patient replica rows
            df_filtered = df_target_ct[gene_level_values.isin(genes_in_data)]

            # 2. Flatten out the gene indices level
            df_reset = df_filtered.reset_index(level="gene")

            # 3. Pivot wider metrics array space into Tidy long layout
            time_cols = [c for c in df_reset.columns if c != "gene"]
            df_long = df_reset.melt(
                id_vars=["gene"],
                value_vars=time_cols,
                var_name="Timepoints",
                value_name="Expression",
            )

            # 4. Generate side-by-side plots layout exactly per notebook spec
            fig, (ax1, ax2) = plt.subplots(
                1, 2,
                figsize=(15, 6),
                gridspec_kw={"width_ratios": [3, 1]},
                constrained_layout=True,
            )

            # --- LEFT PLOT: Expression Lineplot with 95% Confidence Interval ---
            sns.lineplot(
                data=df_long,
                x="Timepoints",
                y="Expression",
                hue="gene",
                marker="o",
                errorbar=("ci", 95), 
                ax=ax1,
            )
            ax1.set_title(f"Pattern {pattern_id}: Expression in {target_ct} Over Time (95% CI)")
            ax1.set_xlabel("Timepoints")
            ax1.set_ylabel("Expression (Log1p)")
            
            col_count = min(6, len(genes_in_data))
            ax1.legend(
                loc="upper center", 
                bbox_to_anchor=(0.5, -0.15), 
                title="Gene",
                ncol=col_count,
                fontsize="small",
                frameon=False 
            )
            ax1.grid(True, linestyle="--", alpha=0.5)

            # --- RIGHT PLOT: Mapper Word Array Shape Curve ---
            try:
                word_arr = mapper[pattern_id].word_array
                x_positions = range(len(word_arr))

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
                ax2.set_xticks(x_positions)
                ax2.set_xticklabels([f"P{k}" for k in x_positions])
                
                y_min, y_max = min(word_arr), max(word_arr)
                if y_min == y_max:
                    ax2.set_ylim(y_min - 1, y_max + 1)
                else:
                    padding = (y_max - y_min) * 0.1
                    ax2.set_ylim(y_min - padding, y_max + padding)
                
                ax2.grid(True, linestyle="--", alpha=0.5)

            except Exception as e:
                ax2.text(
                    0.5, 0.5,
                    f"Error loading mapper[{pattern_id}]:\n{str(e)}",
                    ha="center", va="center",
                )
                ax2.axis("off")

            # --- SAVE THE FIGURE ---
            filename = os.path.join(current_output_dir, f"Subnetwork_{index}_Pattern_{pattern_id}.png")
            plt.savefig(filename, dpi=300, bbox_inches="tight")
            print(f"Saved plot to: {filename}")
            plt.close(fig)  # Closes plot context frame to manage memory bloat in CLI executions

        else:
            print(f"Skipping plot: None of the genes for Subnetwork {index} were found in the expression data.")