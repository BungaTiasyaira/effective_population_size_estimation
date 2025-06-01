# Import modules
import numpy as np
import zarr
import allel
import scipy.signal
from scipy.ndimage import gaussian_filter1d
import dask
import dask.bag as db
from dask.base import compute
from itertools import combinations
from functools import partial
import seaborn as sns
import tskit
import pandas as pd
import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
from collections import Counter, defaultdict
import matplotlib
import matplotlib.pyplot as plt
import csv
import sys
import os
import traceback
import collections
from scoring import get_mapping_of_haplotype_to_mutation, Node, Cluster, Haplotype



tree_path = sys.argv[1]
mutation_position = 35000

def get_paths(tree_path):
    parts = tree_path.split("/")
    seed = parts[0]
    filename = parts[-1].replace(".trees", "")  # Remove extension
    param_parts = filename.split("_")
    allele_freq = param_parts[-1]
    sample_size = "100"
    index = param_parts[1]

    # Construct sample ID path
    sample_ID_path = os.path.join(
        seed,
        "sample_ID",
        allele_freq,
        sample_size,
        f"{filename}_{sample_size}.txt"
    )

    # construct zarr path
    zarr_path = os.path.join(
        seed,
        "ZARR_NEW",
        allele_freq,
        sample_size,
        f"{filename}_{sample_size}.vcz"
    )

    # construct Z_output path
    Z_output = os.path.join(
        seed,
        "PERFECT_TREES",
        allele_freq,
        "Z",
        f"{filename}_{sample_size}.npy"
    )

    # construct dendogram output path
    dendogram_output_path = os.path.join(
        seed,
        "PERFECT_TREES",
        allele_freq,
        "PERFECT_DENDOGRAM",
        f"{filename}_{sample_size}.pdf"
    )

    return sample_ID_path, zarr_path, Z_output, dendogram_output_path, index


sample_ID_path, zarr_path, Z_output, dendogram_output_path, index = get_paths(tree_path)


print(f"Processing tree file: {tree_path}")
print(f"for index: {index}")


def convert(zarr_path, sample_size, mutation_position):

    print("path is", zarr_path)

    # Open as a Dataset
    data = zarr.open_group(zarr_path, mode='r')
    # Stores the ID and genomic position of each variant 
    pos = allel.SortedIndex(data['variant_position'])             
    
    # Extract genotypes for all individuals and convert to haplotypes (sample size x 2)
    gt = data['call_genotype'][:,0:sample_size*2] 
    ht = allel.GenotypeArray(gt).to_haplotypes()    
    
    # Output the frequency of the sweep mutation in the sample
    contains_sweep = pos.locate_range(mutation_position,mutation_position) #finds index of sweep mutation in the array
    sweep = ht[contains_sweep]                           # saves the haplotypes containing the sweep in the variable sweep
    sweep = np.sum(sweep, axis =0)                       #sums up mutation occurences in each haplotypes
    
    samp_freq = np.sum(sweep)/sample_size*2  # finds freq in the entire sample
    

    # This dictionary is used later to color the dendrogram branches according to whether or not the 
    # corresponding sequence contains the sweep mutation
    cols = {}
    for i in range(sample_size*2):
        if sweep[i]:
            cols[i] = "#FF0000" # red 
        else:
            cols[i] = "#808080" # grey 
    
    return cols

cols = convert(zarr_path,100,mutation_position)
dflt_col = "#808080"

# Extract sample node IDs
ts = tskit.load(tree_path)

def read_and_expand_ids(filename):
    expanded_ids = []
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                try:
                    num = int(line)
                    ID1 = 2*num
                    ID2 = ID1 + 1
                    expanded_ids.extend([ID1, ID2])
                except ValueError:
                    print(f"Warning: Skipping non-integer line: '{line}'")
    return expanded_ids

# Example usage
samples = read_and_expand_ids(sample_ID_path)

# Verify they are valid sample nodes
samples = [s for s in samples if ts.node(s).is_sample()]
n = len(samples)

# Compute pairwise distance matrix
dist_matrix = np.zeros((n, n))

def path_length(tree, u, v):
    """Distance from node u up to node v."""
    length = 0
    while u != v:
        parent = tree.parent(u)
        if parent == tskit.NULL:
            break
        length += tree.branch_length(u)
        u = parent
    return length

tree = ts.at(mutation_position)
for i in range(n):
    for j in range(i + 1, n):
        mrca = tree.mrca(samples[i], samples[j])
        dist = path_length(tree, samples[i], mrca) + path_length(tree, samples[j], mrca)
        dist_matrix[i, j] = dist
        dist_matrix[j, i] = dist

# Convert to condensed form and compute linkage
condensed = squareform(dist_matrix)
Z = linkage(condensed, method="average")  


# Ensure all sample indices are in the dict (set to default if missing)
for i in range(n):
    if i not in cols:
        cols[i] = dflt_col

# Compute link colors
link_cols = {}
for i, (left, right) in enumerate(Z[:, :2].astype(int)):
    c1 = link_cols[left] if left >= n else cols[left]
    c2 = link_cols[right] if right >= n else cols[right]
    link_cols[i + n] = c1 if c1 == c2 else dflt_col

# Plot dendrogram
fig, ax = plt.subplots(figsize=(12, 5))
sns.despine(ax=ax, offset=5, bottom=True, top=True)

dd = dendrogram(
    Z,
    link_color_func=lambda x: link_cols[x],
    ax=ax,
    labels=[str(i) for i in range(n)]
)

# Set leaf label colors
for leaf_label in ax.get_xticklabels():
    idx = int(leaf_label.get_text())
    leaf_label.set_color(cols[idx])

plt.tight_layout()


haplotype_to_mutation = get_mapping_of_haplotype_to_mutation(tree_path,
                                                                 sample_ID_path)

mutations = sorted(list(set(haplotype_to_mutation.values())), key=lambda mut: mut.get_generation())
for i, mutation in enumerate(mutations):
    mutation.set_color(matplotlib.colors.to_hex(matplotlib.colormaps['gist_rainbow'](np.linspace(0, 1, len(mutations)))[i]))

haplotypes = [Haplotype(haplotype, haplotype_to_mutation.get(haplotype)) for haplotype in cols]


def add_mutation_info_to_xlabels(plt, haplotype_to_mutation):
    xlabels = plt.gca().get_xticklabels()
    for label in xlabels:
        haplotype = int(label.get_text())
        if haplotype in haplotype_to_mutation:
            mutation = haplotype_to_mutation[haplotype]
            label.set_color(mutation.get_color())
            label.set_text(label.get_text() + " gen " + str(mutation.get_generation()))
    plt.gca().set_xticklabels(xlabels)

fig, ax_dend = plt.subplots(figsize=(30, 8))


# Plot dendrogram and colour branches
fig, ax_dend = plt.subplots(figsize=(20, 8))
ax_dend.set_title('Haplotype clusters',fontsize=24)

dflt_col = "#808080"

node_map = {i: haplotypes[i] for i in range(len(Z)+1)}
clusters = []

for i, row in enumerate(Z[:,:3].astype(int)):
    i += len(Z) + 1
    assert(len(row) == 3)
    l_idx, r_idx, linkage_dist = row
    l, r = node_map[l_idx], node_map[r_idx]

    c1, c2 = [x.get_color() for x in [l, r]]

    this_col = str()
    if c1 == c2:
        this_col = c1
    else:
        this_col = dflt_col

    node = Node(linkage_dist, this_col)
    node.set_left(l)
    node.set_right(r)

    node_map[i] = node

    if ((c1 != dflt_col) and (c2 == dflt_col)):
        clusters.append(Cluster(l))

    if ((c1 == dflt_col) and (c2 != dflt_col)):
        clusters.append(Cluster(r))


sns.despine(ax=ax_dend, offset=5, bottom=True, top=True)
dd = sch.dendrogram(Z,link_color_func=lambda x: node_map[x].get_color(), ax=ax_dend)

ls = []
for leaf, leaf_color in zip(plt.gca().get_xticklabels(), dd["leaves_color_list"]):
    leaf.set_color(haplotypes[int(leaf.get_text())].get_color())
    ls.append(int(leaf.get_text()))

add_mutation_info_to_xlabels(plt, haplotype_to_mutation)

for i in range(len(mutations) + 1):
    ax_dend.plot([]) # hack to allow us to draw legend
ax_dend.set_ylabel('Haplotype age/generations',fontsize=24)
ax_dend.legend([str(int(mutation.get_generation())) for mutation in mutations] + ['mix'], labelcolor=[mutation.get_color() for mutation in mutations] + ['#000000'], handlelength=0, fontsize='xx-small')


os.makedirs(os.path.dirname(dendogram_output_path), exist_ok=True)
os.makedirs(os.path.dirname(Z_output), exist_ok=True)
plt.savefig(dendogram_output_path)
np.save(Z_output, Z)
print(f"saved fig to {dendogram_output_path}")


