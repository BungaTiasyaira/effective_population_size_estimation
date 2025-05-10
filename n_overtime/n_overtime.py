import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from util import parse_log_path
import sys, os

import numpy as np

def plot_nt_with_comparison(stem_filename, N, mu, r, sample_size, color, burnin):
    NUM_SEEDS = 10
    dfs = []
    all_xs = []
    real_ys = []
    artificial_ys = []
    xmin = 0
    xmax = 0

    for i in range(1, NUM_SEEDS + 1):
        # Load CSV
        this_filename = f'{i}/LOG/{i}_{stem_filename}'
        print(this_filename)
        df = pd.read_csv(this_filename)

        df['cycle_adjusted'] = df['cycle'] - 5000

        xmax = max(xmax, df['cycle_adjusted'].iloc[-1])

        dfs.append(df)

    all_xs = np.arange(xmin, xmax + 1)

    for i in range(1, NUM_SEEDS + 1):
        # Load CSV
        df = dfs[i-1]

        # Calculate real nm
        df['nm_real'] = df['Allele Frequency'] * sample_size

        # Calculate real n(t)
        df['n_t_real'] = df['Lineages']

        # Interpolate (spline) the artificial curve
        # Create a finer grid for smoothness

        # Find highest Allele Frequency and compute max_nm
        max_allele_freq = df['Allele Frequency'].max()
        max_nm = round(sample_size * max_allele_freq)

        # Create artificial nm from 1 to max_nm, evenly spaced across same number of points
        nm_artificial = np.linspace(1, max_nm, num=len(df))

        # Calculate artificial n(t)
        n_t_artificial = 4 * N * mu * np.log(1 + (nm_artificial) / (4 * N * mu))

        # plt.plot(df['cycle_adjusted'], df['n_t_real'], marker='o', linestyle='-', label='Real nm')
        # plt.plot(df['cycle_adjusted'], n_t_artificial, marker='x', linestyle='--', color='green', label=f'Artificial nm')

        real_ys.append(gaussian_filter1d(np.interp(all_xs, df['cycle_adjusted'], df['n_t_real']), sigma=2))
        artificial_ys.append(np.interp(all_xs, df['cycle_adjusted'], n_t_artificial))

    avg_n_t_real_smoothed = np.mean(real_ys, axis=0)
    avg_n_t_artificial = np.mean(artificial_ys, axis=0)
    plt.plot(all_xs, avg_n_t_real_smoothed, linestyle='-', color=color, label=f'μ = {mu}; r = {r}')
    plt.plot(all_xs, avg_n_t_artificial, linestyle='--', color=color)
    


if __name__ == "__main__":
    # try:
    #     filename = sys.argv[1] # filename to a log file

    # except IndexError:
    #     print("Usage: python3 analysis.py <filename>")
    #     exit()
    filenames = [
        '1/LOG/1_13_1000_0.00025_2.5e-05_LOG.csv',
        '1/LOG/1_14_1000_0.00025_0.00025_LOG.csv',
        '1/LOG/1_15_1000_0.00025_0.0025_LOG.csv',
        '1/LOG/1_16_1000_0.0025_2.5e-05_LOG.csv',
        '1/LOG/1_17_1000_0.0025_0.00025_LOG.csv'
        ]

    colors = plt.cm.get_cmap('tab10')(np.linspace(0, 1, 10))

    plt.figure(figsize=(10, 6))
    for filename, color in zip(filenames, colors):
        seedless_filename = '_'.join(os.path.basename(filename).split('_')[1:]) # ignore seed value

        # Example usage:
        _, burnin, population_size, mutation_rate, r = parse_log_path(filename)
        plot_nt_with_comparison(seedless_filename, N=population_size*2, mu=mutation_rate, r=r, sample_size=population_size*0.2, color=color, burnin=burnin)

    plt.xlabel('t (Generations)')
    plt.ylabel('Average number of origins ƞ(t)')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("n_overtime/n_overtime.png")  # Saves the plot to a PNG file
    plt.close()  # Closes the plot so it doesn't display in interactive environments