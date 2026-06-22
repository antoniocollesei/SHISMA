import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Load data
data = "results/dataset_execution_times.csv"
df = pd.read_csv(data)

# Parse parameters
df['Genes'] = df['Dataset'].str.extract(r'samples(\d+)').astype(int)
df['Timepoints'] = df['Dataset'].str.extract(r'time(\d+)').astype(int)
df['Patients'] = df['Dataset'].str.extract(r'pats(\d+)').astype(int)

# Ensure Num_Permutations is an integer if it exists
if 'Num_Permutations' in df.columns:
    df['Num_Permutations'] = df['Num_Permutations'].astype(int)
    # FILTER for standard plots to avoid averaging different permutation loads
    df_standard = df[df['Num_Permutations'] == 1000].copy()
else:
    df_standard = df.copy()

sns.set_theme(style="whitegrid")

# ==============================================================================
# PLOT 1: Feature Scalability (Time vs. Genes)
# ==============================================================================
fig, ax = plt.subplots(figsize=(10, 6))
sns.lineplot(
    data=df_standard, 
    x='Genes', 
    y='Elapsed_Time_Seconds', 
    hue='Timepoints',      
    palette='tab10',       # Deep Blue, Orange, Green (High contrast, no yellow)
    marker='o', 
    errorbar=('ci', 95),   
    err_style='bars',      
    err_kws={'capsize': 5, 'elinewidth': 2}, 
    linewidth=2,
    markersize=8,
    ax=ax
)

#plt.title('Feature Scalability (1000 Permutations)', fontsize=14, pad=15)
plt.xlabel('Number of Genes', fontsize=12)
plt.ylabel('Elapsed Time (Seconds, Log Scale)', fontsize=12)
#ax.set_yscale('log')
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
ax.grid(True, which="both", ls="--", linewidth=0.5)
plt.xticks([100, 500, 1000, 5000])

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles, labels=labels, title='Timepoints', title_fontsize=11)
plt.tight_layout()
plt.savefig('results/time_complexity_genes.png', dpi=300)
plt.close()

# ==============================================================================
# PLOT 2: Temporal Scalability (Time vs. Timepoints)
# ==============================================================================
fig, ax = plt.subplots(figsize=(10, 6))
sns.lineplot(
    data=df_standard, 
    x='Timepoints', 
    y='Elapsed_Time_Seconds', 
    hue='Genes',           
    palette='Set1',        # Deep Red, Blue, Green, Purple (Distinct from Plot 1)
    marker='s',            
    errorbar=('ci', 95),   
    err_style='bars',      
    err_kws={'capsize': 5, 'elinewidth': 2},
    linewidth=2,
    markersize=8,
    ax=ax
)

#plt.title('Temporal Scalability (1000 Permutations)', fontsize=14, pad=15)
plt.xlabel('Number of Timepoints', fontsize=12)
plt.ylabel('Elapsed Time (Seconds, Log Scale)', fontsize=12)
ax.set_yscale('log')
ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
plt.xticks([6, 10, 15])    
ax.grid(True, which="both", ls="--", linewidth=0.5)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles, labels=labels, title='Num. Genes', title_fontsize=11)
plt.tight_layout()
plt.savefig('results/time_complexity_timepoints.png', dpi=300)
plt.close()

# ==============================================================================
# PLOT 3: Statistical Scalability (Time vs. Permutations)
# ==============================================================================
if 'Num_Permutations' in df.columns:
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.lineplot(
        data=df, 
        x='Num_Permutations', 
        y='Elapsed_Time_Seconds', 
        hue='Genes', 
        palette='Set1',        # Keeping consistent with Plot 2 since the hue is 'Genes'
        marker='^', 
        errorbar=('ci', 95), 
        err_style='bars', 
        err_kws={'capsize': 5, 'elinewidth': 2},
        linewidth=2,
        markersize=8,
        ax=ax
    )

    #plt.title('Statistical Scalability vs Network Size', fontsize=14, pad=15)
    plt.xlabel('Number of Permutations', fontsize=12)
    plt.ylabel('Elapsed Time (Seconds, Log Scale)', fontsize=12)
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    plt.xticks([10, 100, 500, 1000])
    ax.grid(True, which="both", ls="--", linewidth=0.5)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels, title='Num. Genes', title_fontsize=11)
    plt.tight_layout()
    plt.savefig('results/time_complexity_permutations.png', dpi=300)
    plt.close()