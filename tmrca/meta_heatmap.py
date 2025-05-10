import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore


# Replace 'your_file.csv' with the path to your actual CSV file
csv_file = 'meta.csv'

# Read the CSV file
df = pd.read_csv(csv_file)

# use z scores
df['z_score'] = zscore(df['total_score'])

# get data
heatmap_data = df.pivot_table(
    index='points', 
    columns='threshold', 
    values='z_score'
)

# Plot the heatmap
plt.figure(figsize=(8, 6))
sns.heatmap(heatmap_data, annot=True, cmap="RdYlGn", fmt=".3f")
plt.title("Heatmap of Total Score by Threshold and Points")
plt.xlabel("Peak-Trough Threshold")
plt.ylabel("Smoothing Points")
plt.tight_layout()
plt.savefig("tmrca/heatmap_trial.png")