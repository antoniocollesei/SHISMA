import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

# Assign specific colors for each ppi_type
palette_dict = {
    "quarter": "tab:green",
    "half": "tab:orange",
    "full": "tab:blue"
}

time_data = pd.read_csv("timing_results.csv")

# add column with ratio of time to number of atoms
time_data["rows_cols_ratio"] = time_data["rows"] / time_data["cols"]
time_data["cols_rows_ratio"] = time_data["cols"] / time_data["rows"]

plt.figure(figsize=(10, 6))
sns.lineplot(data=time_data, x="rows", y="runtime_sec", hue="ppi_type", marker="o", palette=palette_dict)
plt.xscale("log")
#plt.yscale("log")
plt.xlabel("Number of Time Series (Log Scale)")
plt.ylabel("Time (seconds)")
plt.legend(title="PPI Size")
plt.grid(True, which="both", ls="--", lw=0.5)
plt.tight_layout()
plt.savefig("plot_rows_vs_runtime.pdf")
plt.savefig("plot_rows_vs_runtime.tiff", dpi=350)

plt.figure(figsize=(10, 6))
sns.lineplot(data=time_data, x="cols", y="runtime_sec", hue="ppi_type", marker="o", palette=palette_dict)
#plt.xscale("log")
#plt.yscale("log")
plt.xlabel("Number of Time Points")
plt.ylabel("Time (seconds)")
plt.title("Timing Results for Different Methods")
plt.legend(title="PPI Size")
plt.grid(True, which="both", ls="--", lw=0.5)
plt.tight_layout()
plt.savefig("plot_cols_vs_runtime.pdf")
plt.savefig("plot_cols_vs_runtime.tiff", dpi=350)

plt.figure(figsize=(10, 6))
sns.lineplot(data=time_data, x="rows_cols_ratio", y="runtime_sec", hue="ppi_type", marker="o", palette=palette_dict)
plt.xscale("log")
#plt.yscale("log")
plt.xlabel("Time Series to Time Points Ratio")
plt.ylabel("Time (seconds)")
plt.title("Timing Results for Different Methods")
plt.legend(title="PPI Size")
plt.grid(True, which="both", ls="--", lw=0.5)
plt.tight_layout()
plt.savefig("plot_rows_cols_ratio_vs_runtime.pdf")
plt.savefig("plot_rows_cols_ratio_vs_runtime.tiff", dpi=350)

time_perm = pd.read_csv("timing_results_permutations.csv")

plt.figure(figsize=(10, 6))
sns.lineplot(data=time_perm, x="perm", y="runtime_sec", hue="ppi_type", marker="o", palette=palette_dict)
#plt.xscale("log")
plt.yscale("log")
plt.xlabel("Number of Permutations")
plt.ylabel("Log Time (seconds)")
plt.legend(title="PPI Size")
plt.grid(True, which="both", ls="--", lw=0.5)
plt.tight_layout()
plt.savefig("plot_perms_vs_runtime.pdf")
plt.savefig("plot_perms_vs_runtime.tiff", dpi=350)