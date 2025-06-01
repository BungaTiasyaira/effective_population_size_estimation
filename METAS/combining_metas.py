import pandas as pd

# Define input files and corresponding method labels
files = {
    'metas/MEOWSS_meta.csv': 'MEOWSS',
    'metas/phasedibd_meta.csv': 'phasedibd',
    'metas/singer_meta.csv': 'SINGER'
}

# Read and label each CSV, then store in a list
dfs = []
for filepath, method_label in files.items():
    df = pd.read_csv(filepath)
    df['method'] = method_label
    dfs.append(df)

# Concatenate all DataFrames
merged_df = pd.concat(dfs, ignore_index=True)

# Save to new CSV
merged_df.to_csv('merged_meta.csv', index=False)
