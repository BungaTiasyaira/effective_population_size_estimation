import pandas as pd

# Load your CSV file (replace 'your_file.csv' with your actual file path)
df = pd.read_csv('metas/merged_meta.csv')

# Create a new column that is the name without the first character (digit)
df['name_modified'] = df['name'].astype(str).str[1:]

# Calculate the difference between n_s and n_p
df['diff_n_s_n_p'] = df['n_s'] - df['n_p']

# Group by 'method' and this new modified name, then get average difference
result = df.groupby(['method', 'name_modified'])['diff_n_s_n_p'].mean().reset_index()

# Rename columns for clarity
result = result.rename(columns={'name_modified': 'name_group', 'diff_n_s_n_p': 'avg_diff_n_s_n_p'})

print(result)

result.to_csv('average_diff_by_method_name.csv', index=False)


