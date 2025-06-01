import pandas as pd
import numpy as np
from decimal import *
import math
from scipy.integrate import simpson
import matplotlib.pyplot as plt

# Likelihood computation function
def p_equation(theta, num_independent_origins, num_mutant_haplotypes):
    with localcontext() as ctx:
        ctx.prec = 32
        theta = Decimal(theta)
        f_x = theta * ctx.ln(1 + Decimal(num_mutant_haplotypes) / theta)
        numerator = ctx.power(f_x, num_independent_origins)
        denominator = math.factorial(num_independent_origins)
        e_term = ctx.power(Decimal(math.e), -f_x)
        return float((numerator / denominator) * e_term)

def plot_likelihood_curve(ax, df, color, label_prefix, chosen_n, mu):
    theta_values = np.linspace(0.01, 300, 10000)
    Ne_values = theta_values / (4 * mu)
    pdf_matrix = []

    for _, row in df.iterrows():
        pdf = np.array([p_equation(theta, int(row[chosen_n]), int(row['k'])) for theta in theta_values])
        normalized_pdf = pdf / simpson(pdf, theta_values)
        pdf_matrix.append(normalized_pdf)

    combined_pdf = np.prod(pdf_matrix, axis=0)
    combined_pdf /= simpson(combined_pdf, theta_values)
    combined_pdf /= np.max(combined_pdf)

    threshold = np.exp(-3.84 / 2)
    inside_CI = theta_values[combined_pdf >= threshold]

    theta_lower = inside_CI[0]
    theta_upper = inside_CI[-1]
    theta_mle = theta_values[np.argmax(combined_pdf)]

    Ne_lower = theta_lower / (4 * mu)
    Ne_upper = theta_upper / (4 * mu)
    Ne_mle = theta_mle / (4 * mu)

    ax.plot(Ne_values, combined_pdf, label=f"{label_prefix}", color=color, linewidth=2)
    ax.axvline(Ne_lower, color=color, linestyle='--', linewidth=1)
    ax.axvline(Ne_upper, color=color, linestyle='--', linewidth=1)

    print(f"{label_prefix} Combined 95% CI (theta): [{theta_lower:.3f}, {theta_upper:.3f}]")
    print(f"{label_prefix} MLE (theta): {theta_mle:.3f}")

# Load data
df_input = pd.read_csv("tmrca/get_likelihood_input.csv")
fig1, ax1 = plt.subplots(figsize=(12, 7))
fig2, ax2 = plt.subplots(figsize=(12, 7))

param_dict = {
    "param1": ["_13_1000_0.00025_2.5e-05_0.5_100", "\u03c1 = 0.1, t = 4, AF = 0.5", "n4", '#669fff'],
    "param2": ["_14_1000_0.00025_0.00025_0.5_100", "\u03c1 = 1.0, t = 4, AF = 0.5", "n4", '#005fff'],
    "param3": ["_15_1000_0.00025_0.0025_0.5_100", "\u03c1 = 10, t = 4, AF = 0.5", "n4", '#003080'],
    "param4": ["_13_1000_0.00025_2.5e-05_0.7_100", "\u03c1 = 0.1, t = 4, AF = 0.7", "n4", '#c466ff'],
    "param5": ["_14_1000_0.00025_0.00025_0.7_100", "\u03c1 = 1.0, t = 4, AF = 0.7", "n4", '#9c00ff'],
    "param6": ["_15_1000_0.00025_0.0025_0.7_100", "\u03c1 = 10, t = 4, AF = 0.7", "n4", '#5e0099'],
    "param7": ["_16_1000_0.0025_2.5e-05_0.5_100", "\u03c1 = 0.1, t = 1, AF = 0.5", "n1", '#669fff'],
    "param8": ["_17_1000_0.0025_0.00025_0.5_100", "\u03c1 = 1.0, t = 2, AF = 0.5", "n2", '#003080'],
    "param9": ["_16_1000_0.0025_2.5e-05_0.7_100", "\u03c1 = 0.1, t = 1, AF = 0.7", "n1", '#c466ff'],
    "param10": ["_17_1000_0.0025_0.00025_0.7_100", "\u03c1 = 1.0, t = 2, AF = 0.7", "n2", '#5e0099']
}

df = {}
for i in range(1, 11):
    key = f"param{i}"
    suffix, label, chosen_n, color = param_dict[key]
    mu = 0.00025 if i <= 6 else 0.0025
    ax = ax1 if i <= 6 else ax2

    df[i] = df_input[df_input["label"].str.endswith(suffix)].copy()
    print(suffix)
    plot_likelihood_curve(ax, df[i], color=color, label_prefix=label, chosen_n=chosen_n, mu=mu)

for ax in [ax1, ax2]:
    ax.set_xlabel("Effective population size (Ne)")
    ax.set_ylabel("Normalized Likelihood")
    ax.set_title("Likelihood Curves with 95% Confidence Intervals")
    ax.set_xlim(0, 10000)
    ax.grid(True)
    ax.legend(loc='upper right', fontsize=14)

plt.tight_layout()
fig1.savefig("tmrca/MEOWSS_0.png")
fig2.savefig("tmrca/MEOWSS_1.png")
