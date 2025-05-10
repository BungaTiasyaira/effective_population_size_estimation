import pandas as pd
import matplotlib.pyplot as plt

# Replace 'your_file.csv' with the path to your actual CSV file
csv_file = 'meta.csv'

# Read the CSV file
df = pd.read_csv(csv_file)

# Plotting
plt.figure(figsize=(10, 6))

# Group by 'points' and calculate the mean of total_score for each point value
for allele_freq in ["_0.5_", "_0.7_"]:
    vcf_points = df[df['vcf'].str.contains(allele_freq)]
    print(vcf_points)
    means = vcf_points.groupby('points')['total_score'].mean().reset_index()

    for threshold in vcf_points['threshold'].unique():
        threshold_points = vcf_points.loc[vcf_points['threshold'] == threshold]
        means = threshold_points.groupby('points')['total_score'].mean().reset_index()
        print(means)
        plt.plot(means['points'], means['total_score'], marker='o', linestyle='-', label=f'graph {allele_freq} {threshold}')

    # means = vcf_points.groupby('threshold')['total_score'].mean().reset_index()
    # print(means)
    # plt.plot(means['threshold'], means['total_score'], marker='o', linestyle='-', label=f'graph {allele_freq}')

plt.xlabel('Points')
plt.ylabel('Average Total Score')
plt.title('Average Total Score vs Points')
plt.grid(True)
plt.tight_layout()
plt.legend()

# Show the plot
plt.savefig('meta2.png')