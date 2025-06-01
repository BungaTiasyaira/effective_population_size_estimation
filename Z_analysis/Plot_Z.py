import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Define method colors
method_colors = {
    "MEOWSS": "#f99bd2",
    "SINGER": "#9bc3f9",
    "phasedibd": "#9bf9b2"
}

x_labels = {
    13: "θ = 1.0; ρ = 0.1",
    14: "θ = 1.0; ρ = 1.0",
    15: "θ = 1.0; ρ = 10",
    16: "θ = 10; ρ = 0.1",
    17: "θ = 10; ρ = 1.0"
}

# Load and preprocess data
df = pd.read_csv("Z_analysis/CophCorr.csv")  # Update path
melted = df.melt(id_vars=["seed", "simulation", "AF"], 
                 value_vars=["SINGER", "phasedibd", "MEOWSS"], 
                 var_name="method", value_name="score")
melted = melted.dropna(subset=["score"])

# Aggregate mean and std
summary = (melted.groupby(["AF", "simulation", "method"])
                 .agg(mean_score=("score", "mean"), std_score=("score", "std"))
                 .reset_index())

# Plot for each allele frequency
for af, group in summary.groupby("AF"):
    simulations = sorted(group["simulation"].unique())
    methods = ["SINGER", "phasedibd", "MEOWSS"]

    bar_width = 0.25
    x = np.arange(len(simulations))  # one group per simulation

    plt.figure(figsize=(10, 6))

    for i, method in enumerate(methods):
        data = group[group["method"] == method]
        means = []
        stds = []
        for sim in simulations:
            row = data[data["simulation"] == sim]
            if not row.empty:
                means.append(row["mean_score"].values[0])
                stds.append(row["std_score"].values[0])
            else:
                means.append(np.nan)
                stds.append(0)

        x_pos = x + (i - 1) * bar_width
        plt.bar(x_pos, means, yerr=stds, width=bar_width,
                label=method, color=method_colors.get(method, "gray"), capsize=5)

    plt.xticks(x, [f"{x_labels.get(s)}" for s in simulations], fontsize=12)
    plt.yticks(fontsize=14)
    plt.title(f"Clustering Performance (allele freq. = {af})")
    plt.ylabel("Mean Cophenetic Correlation Coefficient", fontsize=12)
    plt.legend(title="Method")
    plt.tight_layout()
    plt.savefig(f"Z_analysis/CophCorr_{af}.png")
    plt.close()

from scipy.stats import f_oneway, ttest_ind
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import warnings
warnings.filterwarnings("ignore")

# Melted data: one row per seed per method
results = []

for (af, sim), group in melted.groupby(["AF", "simulation"]):
    # Filter out methods with no data
    available_methods = group["method"].unique()
    scores_by_method = [group[group["method"] == m]["score"].values for m in available_methods]

    if len(scores_by_method) >= 2:
        # ANOVA
        fval, pval = f_oneway(*scores_by_method)
        results.append((af, sim, "ANOVA", "all", pval))

        # Pairwise t-tests
        for i in range(len(available_methods)):
            for j in range(i + 1, len(available_methods)):
                m1 = available_methods[i]
                m2 = available_methods[j]
                s1 = group[group["method"] == m1]["score"]
                s2 = group[group["method"] == m2]["score"]
                if len(s1) > 1 and len(s2) > 1:
                    _, pair_pval = ttest_ind(s1, s2, equal_var=False)
                    results.append((af, sim, "t-test", f"{m1} vs {m2}", pair_pval))

# Convert to DataFrame for inspection
stat_df = pd.DataFrame(results, columns=["AF", "Simulation", "Test", "Comparison", "p-value"])

# Optional: apply multiple testing correction
from statsmodels.stats.multitest import multipletests
stat_df["adjusted_p"] = multipletests(stat_df["p-value"], method="fdr_bh")[1]

# Filter significant results (e.g., p < 0.05)
significant = stat_df[stat_df["adjusted_p"] < 0.05]
print(significant)

# Add significance flag
stat_df["significant"] = stat_df["adjusted_p"] < 0.05

# Save to CSV
stat_df.to_csv("Z_analysis/statistical_significance_results.csv", index=False)