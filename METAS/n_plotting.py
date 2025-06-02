import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

method_color_map = {
    "phasedibd": {
        "n_d_1": "#afeccc", "n_d_2": "#79cea1", "n_d_3": "#47a673", "n_d_4": "#1e7345",
        "Ne_d_1": "#afeccc", "Ne_d_2": "#79cea1", "Ne_d_3": "#47a673", "Ne_d_4": "#1e7345"
    },
    "SINGER": {
        "n_d_1": "#6f2045", "n_d_2": "#668ec5", "n_d_3": "#3a649e", "n_d_4": "#113567",
        "Ne_d_1": "#6f2045", "Ne_d_2": "#668ec5", "Ne_d_3": "#3a649e", "Ne_d_4": "#113567"
    },
    "MEOWSS": {
        "n_d_1": "#f1c2d8", "n_d_2": "#d692b2", "n_d_3": "#a44f77", "n_d_4": "#6f2045",
        "Ne_d_1": "#f1c2d8", "Ne_d_2": "#d692b2", "Ne_d_3": "#a44f77", "Ne_d_4": "#6f2045"
    }
}

ground_truth_color_map = {
    "n_s": "#f9a323", "n_p": "#ef383a",
    "Ne_s": "#f9a323", "Ne_p": "#ef383a"
}

legend_label_map = {
    "MEOWSS_n_d_1": "MEOWSS; t = 1",
    "MEOWSS_n_d_2": "MEOWSS; t = 2",
    "MEOWSS_n_d_3": "MEOWSS; t = 3",
    "MEOWSS_n_d_4": "MEOWSS; t = 4",
    "phasedibd_n_d_1": "PhasedIBD; t = 1",
    "phasedibd_n_d_2": "PhasedIBD; t = 2",
    "phasedibd_n_d_3": "PhasedIBD; t = 3",
    "phasedibd_n_d_4": "PhasedIBD; t = 4",
    "SINGER_n_d_1": "SINGER; t = 1",
    "SINGER_n_d_2": "SINGER; t = 2",
    "SINGER_n_d_3": "SINGER; t = 3",
    "SINGER_n_d_4": "SINGER; t = 4",
    "GT_n_s": "GT-Sample",
    "GT_n_p": "GT-Population",
    "MEOWSS_Ne_d_1": "MEOWSS; t = 1",
    "MEOWSS_Ne_d_2": "MEOWSS; t = 2",
    "MEOWSS_Ne_d_3": "MEOWSS; t = 3",
    "MEOWSS_Ne_d_4": "MEOWSS; t = 4",
    "phasedibd_Ne_d_1": "PhasedIBD; t = 1",
    "phasedibd_Ne_d_2": "PhasedIBD; t = 2",
    "phasedibd_Ne_d_3": "PhasedIBD; t = 3",
    "phasedibd_Ne_d_4": "PhasedIBD; t = 4",
    "SINGER_Ne_d_1": "SINGER; t = 1",
    "SINGER_Ne_d_2": "SINGER; t = 2",
    "SINGER_Ne_d_3": "SINGER; t = 3",
    "SINGER_Ne_d_4": "SINGER; t = 4",
    "GT_Ne_s": "GT-Sample",
    "GT_Ne_p": "GT-Population"
}

# Load CSV
df = pd.read_csv("metas/merged_meta.csv")

# Convert relevant columns to numeric
non_numeric = ['name', 'method']
numeric_cols = [col for col in df.columns if col not in non_numeric]
df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce')

# Parse 'name' to extract sim_id and allele_freq
df[['seed', 'sim_id', 'pop_size', 'mut_rate', 'rec_rate', 'allele_freq', 'something']] = df['name'].str.split('_', expand=True)
df['sim_id'] = df['sim_id'].astype(int)
df['allele_freq'] = df['allele_freq'].astype(float)

# Group by sim_id and allele_freq, average numeric columns
grouped_mean = df.groupby(['method', 'sim_id', 'allele_freq'])[numeric_cols].mean()
grouped_std = df.groupby(['method', 'sim_id', 'allele_freq'])[numeric_cols].std()

grouped_mean = grouped_mean.reset_index()
grouped_std = grouped_std.reset_index()

# Define the full sets of columns
full_n_cols = ['n_d_1','n_d_2','n_d_3','n_d_4','n_s','n_p']
full_Ne_cols = ['Ne_d_1','Ne_d_2','Ne_d_3','Ne_d_4','Ne_s','Ne_p']

# Create identifier for x-axis
grouped_mean['sim_label'] = grouped_mean['sim_id'].astype(str) + "_" + grouped_mean['allele_freq'].astype(str)

# Define your simulation ID groups, columns to plot, and colors for each column
plot_configs = [
    {
        'sim_ids': [13, 14, 15],
        'n_cols': ['n_d_3','n_d_4', 'n_s','n_p'],
        'Ne_cols': ['Ne_d_3','Ne_d_4', 'Ne_s','Ne_p'],
        'prefix': 'sim_13_14_15'
    },
    {
        'sim_ids': [16, 17],
        'n_cols': full_n_cols,
        'Ne_cols': full_Ne_cols,
        'prefix': 'sim_16_17'
    }
]

def plot_using_matplotlib(ax, data, data_std, colors):
    x = np.arange(len(data.columns))  # positions for each group on x-axis
    n_bars = len(data.index)          # number of bars in each group
    total_width = 0.8                        # total width reserved per group
    bar_width = total_width / n_bars        # width of each individual bar

    for i, (label, row) in enumerate(data.iterrows()):
        offset = (i - n_bars / 2) * bar_width + bar_width / 2
        std_row = data_std.loc[label]

        ax.bar(x + offset, row.values, width=bar_width, color=colors[i], label=label, yerr=std_row.values,
        capsize=5)

    ax.set_xticks(x)
    ax.set_xticklabels(data.columns, rotation=45, ha='right')



for config in plot_configs:
    sim_group_mean = grouped_mean[grouped_mean['sim_id'].isin(config['sim_ids'])]
    sim_group_std = grouped_std[grouped_std['sim_id'].isin(config['sim_ids'])]

    # === For n plots (per AF) ===
    for af in sorted(sim_group_mean['allele_freq'].unique()):
        sub_mean = sim_group_mean[sim_group_mean['allele_freq'] == af].copy()
        sub_std = sim_group_std[sim_group_std['allele_freq'] == af].copy()

        bar_data_n, bar_data_n_std = pd.DataFrame(), pd.DataFrame()
        xticks_labels_n = []
        colors_n = []

        for method in sub_mean['method'].unique():
            mean_group = sub_mean[sub_mean['method'] == method]
            std_group = sub_std[sub_std['method'] == method]

            for _, row in mean_group.iterrows():
                label = f"{row['sim_id']}_{af}"
                std_row = std_group[(std_group['sim_id'] == row['sim_id']) & (std_group['allele_freq'] == row['allele_freq'])].squeeze()
                xticks_labels_n.append(label)

                # Collect n-values per method
                for col in config['n_cols']:
                    if 'd_' in col:
                        bar_data_n.loc[f"{method}_{col}", label] = row[col]
                        bar_data_n_std.loc[f"{method}_{col}", label] = std_row[col]
                        colors_n.append(method_color_map[method][col])
                    elif 's' in col or 'p' in col:
                        bar_data_n.loc[f"GT_{col}", label] = row[col]
                        bar_data_n_std.loc[f"GT_{col}", label] = std_row[col]
                        colors_n.append(ground_truth_color_map[col])

        # Reorder bar_data_n rows
        def put_in_desired_order(data):
            desired_order = []
            desired_order += [idx for idx in data.index if idx.startswith('MEOWSS_')]
            desired_order += [idx for idx in data.index if idx.startswith('phasedibd_')]
            desired_order += [idx for idx in data.index if idx.startswith('SINGER_')]
            desired_order += [idx for idx in data.index if idx.startswith('GT_')]
            return data.loc[desired_order]

        bar_data_n = put_in_desired_order(bar_data_n)
        bar_data_n_std = put_in_desired_order(bar_data_n_std)
        
        colors_n = [
            method_color_map.get(row.split('_')[0], ground_truth_color_map).get('_'.join(row.split('_')[1:]), '#000000')
            if not row.startswith('GT') else ground_truth_color_map.get(row.replace("GT_", ""), '#000000')
            for row in bar_data_n.index
        ]

        bar_data_n.rename(index=legend_label_map, inplace=True)
        bar_data_n_std.rename(index=legend_label_map, inplace=True)



        # === Plot n values ===
        fig, ax = plt.subplots(figsize=(20, 12))
        plot_using_matplotlib(ax, bar_data_n, bar_data_n_std, colors_n)
        # ax = bar_data_n.T.plot(kind='bar', yerr=bar_data_n_std.values, figsize=(12,6), color=colors_n, title=f"Number of Origins at AF={af}")
        added_lines = set()

        for sim_id in sub_mean['sim_id'].unique():
            if sim_id in [13, 14, 15]:
                if af == 0.5 and 'theta1_af0.5' not in added_lines:
                    ax.axhline(y=4.615, color='red', linestyle='--', linewidth=1.5, label='n modelled (θ=1, AF=0.5)')
                    added_lines.add('theta1_af0.5')
                elif af == 0.7 and 'theta1_af0.7' not in added_lines:
                    ax.axhline(y=4.95, color='red', linestyle='--', linewidth=1.5, label='n modelled (θ=1, AF=0.7)')
                    added_lines.add('theta1_af0.7')
            elif sim_id in [16, 17]:
                if af == 0.5 and 'theta10_af0.5' not in added_lines:
                    ax.axhline(y=23.13, color='red', linestyle='--', linewidth=1.5, label='n modelled (θ=10, AF=0.5)')
                    added_lines.add('theta10_af0.5')
                elif af == 0.7 and 'theta10_af0.7' not in added_lines:
                    ax.axhline(y=26.46, color='red', linestyle='--', linewidth=1.5, label='n modelled (θ=10, AF=0.7)')
                    added_lines.add('theta10_af0.7')

        plt.ylabel("Independent Origins")
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(f"metas/graphs/combined_n_{config['prefix']}_af_{af}.png")
        plt.close()

    # === For Ne plot (combine all AFs into one) ===
    bar_data_Ne, bar_data_Ne_std = pd.DataFrame(), pd.DataFrame()
    xticks_labels_Ne = []
    colors_Ne = []

    for af in sorted(sim_group_mean['allele_freq'].unique()):
        sub_mean = sim_group_mean[sim_group_mean['allele_freq'] == af].copy()
        sub_std = sim_group_std[sim_group_std['allele_freq'] == af].copy()

        for method in sub_mean['method'].unique():
            mean_group = sub_mean[sub_mean['method'] == method]
            std_group = sub_std[sub_std['method'] == method]

            for _, row in mean_group.iterrows():
                label = f"{row['sim_id']}_{af}"
                std_row = std_group[(std_group['sim_id'] == row['sim_id']) & (std_group['allele_freq'] == row['allele_freq'])].squeeze()
                xticks_labels_Ne.append(label)

                # Collect Ne-values per method
                for col in config['Ne_cols']:
                    if 'd_' in col:
                        bar_data_Ne.loc[f"{method}_{col}", label] = row[col]
                        bar_data_Ne_std.loc[f"{method}_{col}", label] = std_row[col]
                        colors_Ne.append(method_color_map[method][col])
                    elif 's' in col or 'p' in col:
                        bar_data_Ne.loc[f"GT_{col}", label] = row[col]
                        bar_data_Ne_std.loc[f"GT_{col}", label] = std_row[col]
                        colors_Ne.append(ground_truth_color_map[col])

    # Reorder bar_data_Ne rows
    bar_data_Ne = put_in_desired_order(bar_data_Ne)
    bar_data_Ne_std = put_in_desired_order(bar_data_Ne_std)

    colors_Ne = [
        method_color_map.get(row.split('_')[0], ground_truth_color_map).get('_'.join(row.split('_')[1:]), '#000000')
        if not row.startswith('GT') else ground_truth_color_map.get(row.replace("GT_", ""), '#000000')
        for row in bar_data_Ne.index
    ]

    bar_data_Ne.rename(index=legend_label_map, inplace=True)
    bar_data_Ne_std.rename(index=legend_label_map, inplace=True)


    # === Plot Ne values ===
    plt.rcParams.update({'font.size': 20})
    fig, ax = plt.subplots(figsize=(20, 12))

    plot_using_matplotlib(ax, bar_data_Ne, bar_data_Ne_std, colors_Ne)

    # Add horizontal line for true Ne
    ax.axhline(y=1000, color='red', linestyle='--', linewidth=1.5, label='True Ne = 1000')

    # Format ticks, labels, legend
    ax.set_ylabel("Effective Population Size")
    ax.legend(fontsize=16, loc='upper left')
    plt.tight_layout()
    plt.savefig(f"metas/graphs/combined_Ne_{config['prefix']}.png")
    plt.close()



