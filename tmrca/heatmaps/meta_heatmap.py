import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore

# Read the CSV file
csv_file = 'meta.csv'
df = pd.read_csv(csv_file)

# Extract vcf_group (1st underscore-separated group after prefix)
df['vcf_group'] = df['vcf'].str.extract(r'^[^_]+_([^_]+)')

# Extract allele frequency (2nd-to-last underscore-separated part)
df['allele_freq'] = df['vcf'].str.extract(r'(?:[^_]+_){5}([^_]+)')

# Optionally convert to float for numeric sorting/labeling
df['allele_freq'] = df['allele_freq'].astype(float)

# Use z-scores within each (vcf_group, allele_freq)
df['z_score'] = df.groupby(['vcf_group', 'allele_freq'])['total_score'].transform(zscore)

# Group by both
grouped = df.groupby(['vcf_group', 'allele_freq'])

# Create a heatmap for each subgroup
for (group_name, allele_freq), group in grouped:

    heatmap_data = group.pivot_table(
        index='points', 
        columns='threshold', 
        values='z_score'
    )

    plt.figure(figsize=(8, 6))
    sns.heatmap(heatmap_data, annot=True, cmap="RdYlGn", fmt=".3f")
    plt.title(f"Z-Score Heatmap\nGroup: {group_name} | Allele Freq: {allele_freq}")
    plt.xlabel("Peak-Trough Threshold")
    plt.ylabel("Smoothing Points")
    plt.tight_layout()
    # Use underscores instead of dots in file name to avoid filesystem issues
    safe_freq = str(allele_freq).replace('.', '_')
    plt.savefig(f"tmrca/heatmap_{group_name}_{safe_freq}.png")
    plt.close()

