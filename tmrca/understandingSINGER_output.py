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
from analysis import run_tmrca_analysis, get_paths, get_vcf_basename

# === Load your tree sequence ===
tree_path = sys.argv[1]
index = int(sys.argv[2])

print(f"Processing tree file: {tree_path}")
print(f"for index: {index}")

def get_zarr_from_trees (tree_path, index):
    parts = tree_path.split("/")
    filename = parts[1]
    params = filename.split("_")
    seed = params[0]
    sample_size = params[-1]
    AF = params[-2]

    zarr_path = os.path.join(
        seed,
        "ZARR_NEW",
        AF,
        sample_size,
        f"{filename}.vcz"
    )

    Z_output = os.path.join(
        seed,
        "SINGER",
        AF,
        "Z",
        f"{filename}_{index}.npy"
    )
    os.makedirs(os.path.dirname(Z_output), exist_ok=True)

    dendogram_output = os.path.join(
        seed,
        "SINGER",
        AF,
        "Dendograms",
        f"{filename}_{index}.pdf"
    )
    os.makedirs(os.path.dirname(dendogram_output), exist_ok=True)

    return zarr_path, Z_output, dendogram_output

zarr_path, Z_output, dendogram_output = get_zarr_from_trees(tree_path, index)


_, _, _, trees_path, sample_ID_path, population_size, mutation_rate = get_paths(zarr_path)


mutation_position = 35001

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
samples = ts.samples()
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
np.save(Z_output, Z)

# Ensure all sample indices are in the dict (set to default if missing)
for i in range(n):
    if i not in cols:
        cols[i] = dflt_col

run_tmrca_analysis(Z, cols, trees_path, sample_ID_path, dendogram_output, population_size, mutation_rate, "singer_meta.csv", get_vcf_basename(zarr_path))
