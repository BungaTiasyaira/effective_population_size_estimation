import phasedibd as ibd
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
import scipy.cluster.hierarchy as sch
from collections import Counter, defaultdict
import matplotlib
import csv
import sys
import os
import traceback
import collections
import pandas as pd
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
from analysis import run_tmrca_analysis, get_vcf_basename
import zarr


def get_paths(vcf_path):
    parts = vcf_path.split("/")
    seed = parts[0]
    allele_freq = parts[2]
    sample_size = parts[3]
    filename = parts[-1].replace(".vcf", "")
    filename_parts = filename.split("_")
    index = filename_parts[1]
    mutation_rate = filename_parts[3]
    rec_rate = filename_parts[4]
    for_trees = filename[:-4]

    # get trees path 
    trees_path = os.path.join(
        seed,
        "TREES",
        f"{for_trees}.trees"
    )

    # get sample id path
    sample_ID_path = os.path.join(
        seed,
        "sample_ID",
        allele_freq,
        sample_size,
        f"{filename}.txt"
    )

    # get zarr_path 
    zarr_path = os.path.join(
        seed,
        "ZARR_NEW",
        allele_freq,
        sample_size,
        f"{filename}.vcz"
    )

    # get dendogram_path
    dendogram_path = os.path.join(
        seed,
        "PHASEIBD",
        allele_freq,
        "dendograms",
        f"{filename}.pdf"
    )
    os.makedirs(os.path.dirname(dendogram_path), exist_ok=True)

    # get Z_output
    Z_output = os.path.join(
        seed,
        "PHASEIBD",
        allele_freq,
        "Z",
        f"{filename}.npy"
    )
    os.makedirs(os.path.dirname(Z_output), exist_ok=True)

    # get csv_output path
    csv_output = os.path.join(
        seed,
        "PHASEIBD",
        allele_freq,
        "csv_output",
        f"{filename}.csv"
    )
    os.makedirs(os.path.dirname(csv_output), exist_ok=True)

    # get csv_unique path
    csv_unique = os.path.join(
        seed,
        "PHASEIBD",
        allele_freq,
        "csv_unique",
        f"{filename}.csv"
    )
    os.makedirs(os.path.dirname(csv_unique), exist_ok=True)

    return trees_path, sample_ID_path, zarr_path, dendogram_path, Z_output, mutation_rate, csv_output, csv_unique, rec_rate, index


# uses vcz to get cols (for dendogram) and sweep
def convert(zarr_path, sample_size, mutation_position):
    """
    extracts haplotypes from zarr file and identifies which haplotype has the sweep mutation 

    output: cols (for tmrca analysis), sweep (for phasedibd function)
    """

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
    
    return cols, sweep, ht, pos


def count_haplotype_differences(ht, zarr_path, start_bp, end_bp, hap1_idx, hap2_idx, pos):
    """
    Count allele differences between two haplotypes within a genomic region.

    output: differences in SNPs between a pair of haplotypes 

    """
   
    try:
        region_slice = pos.locate_range(start_bp, end_bp)
        htx = ht[region_slice]

        if htx.shape[0] == 0:
            return 0

        hap1 = htx[:, hap1_idx]
        hap2 = htx[:, hap2_idx]
        diff = np.sum(hap1 != hap2)

        print(f"Comparing haplotypes {hap1_idx} vs {hap2_idx} in region {start_bp}-{end_bp}, variants: {htx.shape[0]}, differences: {diff}")
    except Exception:
        diff = 0
        print("Error in count_haplotype_differences (returning 0):")
        print(traceback.format_exc())

    return diff

# function to perform ibd analysis using phaseibd
def run_phasedIBD(zarr_path, vcf_path, csv_output, csv_unique, ht, mutation_rate, rec_rate, sweep, pos, index):
    """
    Function that performs hierarchichal clustering using identified IBD segments within a region around the focal mutation 

    output: Z distance linkage matrix

    """
    # extract haplotypes
    haplotypes = ibd.VcfHaplotypeAlignment(vcf_path)

    # perform the IBD analysis
    tpbwt = ibd.TPBWTAnalysis()

    # output is a pandas DataFrame with all the IBD segments
    ibd_results = tpbwt.compute_ibd(haplotypes, L_m = 10, L_f = 0.0001, use_phase_correction = False)

    # Compute hap1_name and hap2_name
    ibd_results['hap1_name'] = ibd_results['id1'] * 2 + ibd_results['id1_haplotype']
    ibd_results['hap2_name'] = ibd_results['id2'] * 2 + ibd_results['id2_haplotype']


    # filter and clip regions according to index
    if int(index) == 13:
        region_start = 32300
        region_end = 37700
        tmrca_threshold = 800
        sigma = 300
    elif int(index) == 14:
        region_start = 33000
        region_end = 37000
        tmrca_threshold = 300
        sigma = 300
    elif int(index) == 15:
        region_start = 34000
        region_end = 36000
        tmrca_threshold = 80
        sigma = 300
    elif int(index) == 16:
        region_start = 32300
        region_end = 37700
        tmrca_threshold = 100
        sigma = 600
    elif int(index) == 17:
        region_start = 32300
        region_end = 37700
        tmrca_threshold = 150
        sigma = 600

    # parameters that stay the same
    bonus = 0.3  
    penalty = 1.7 
    scaling_factor = 1000

    ibd_filtered = ibd_results[(ibd_results['end_bp'] > region_start) & (ibd_results['start_bp'] < region_end)].copy()
    ibd_filtered['start_bp'] = ibd_filtered['start_bp'].clip(lower=region_start)
    ibd_filtered['end_bp'] = ibd_filtered['end_bp'].clip(upper=region_end)

    # save output to CSV
    ibd_filtered.to_csv(csv_output, index=False)

    # create new column in csv for the differences within each region
    ibd_filtered['differences'] = ibd_filtered.apply(
        lambda row: count_haplotype_differences(
            ht,
            zarr_path,
            int(row['start_bp']),
            int(row['end_bp']),
            int(row['hap1_name']),
            int(row['hap2_name']),
            pos
        ),
        axis=1
    )

    # Make unordered pairs: ensure each (a, b) is always (min, max)
    observed_pairs = set(tuple(sorted((a, b))) for a, b in zip(ibd_filtered['hap1_name'], ibd_filtered['hap2_name']))

    # Generate all possible unique unordered haplotype pairs (excluding self-pairs)
    all_haplotypes = list(range(200))
    all_possible_pairs = set(combinations(all_haplotypes, 2))
    # Compare sets, print to see if any missing pairs
    existing_count = len(observed_pairs)
    missing_count = len(all_possible_pairs - observed_pairs)
    print(f"Existing pairs: {existing_count}")
    print(f"Missing pairs: {missing_count}")

    # Always sort haplotype pairs to treat them as unordered
    ibd_filtered['pair'] = ibd_filtered.apply(lambda row: tuple(sorted((row['hap1_name'], row['hap2_name']))), axis=1)

    # Compute segment length
    ibd_filtered['segment_length'] = ibd_filtered['end_bp'] - ibd_filtered['start_bp']

    # Define weighting function (Gaussian around focal mutation)
    def ibd_weight(start, end):
        """
        Makes it so that ibd segments closer to the focal mutation have a greater weight than those further away

        """
        focal = 35001

        center = (start + end) / 2
        dist = abs(center - focal)
        return np.exp(- (dist**2) / (2 * sigma**2))
    
    # Apply weights to each IBD segment
    ibd_filtered['weight'] = ibd_filtered.apply(lambda row: ibd_weight(row['start_bp'], row['end_bp']), axis=1)
    ibd_filtered['weighted_length'] = ibd_filtered['segment_length'] * ibd_filtered['weight']
    ibd_filtered['weighted_diff'] = ibd_filtered['differences'] * ibd_filtered['weight']

    # Aggregate using weighted values
    ibd_sums = ibd_filtered.groupby('pair')[['weighted_length', 'weighted_diff']].sum().reset_index()

    # Split the 'pair' into separate haplotype columns
    ibd_sums[['hap1_name', 'hap2_name']] = pd.DataFrame(ibd_sums['pair'].tolist(), index=ibd_sums.index)

    # Clean up
    ibd_sums = ibd_sums.drop(columns='pair').rename(columns={
        'weighted_length': 'IBD_length',
        'weighted_diff': 'total_differences'
    })

    # Get all 19,900 possible pairs, put in all_pairs_df
    all_pairs = list(combinations(range(200), 2))
    all_pairs_df = pd.DataFrame(all_pairs, columns=['hap1_name', 'hap2_name'])

    # Merge, filling missing pairs with IBD_length = 0
    final_df = pd.merge(all_pairs_df, ibd_sums, on=['hap1_name', 'hap2_name'], how='left').fillna({'IBD_length': 100})
    final_df['IBD_length'] = final_df['IBD_length'].astype(int)
    final_df.to_csv(csv_unique, index=False)

    # Get sorted list of all haplotype names
    haplotypes = sorted(set(final_df['hap1_name']).union(final_df['hap2_name']))
    hap_index = {hap: i for i, hap in enumerate(haplotypes)}
    n = len(haplotypes)
    
    # --- mutation-aware distance computation with smooth weighting ---

    # Initialize TMRCA-based distance matrix
    distance_matrix = np.zeros((n, n))
    np.fill_diagonal(distance_matrix, 0)


    # Ensure mutation_status is computed before this block
    mutation_status = {hap: int(val > 0) for hap, val in enumerate(sweep)}

    # Compute TMRCA-based distances
    for _, row in final_df.iterrows():
        i = hap_index[row['hap1_name']]
        j = hap_index[row['hap2_name']]
        ibd_len = float(row['IBD_length'])
        k = float(row['total_differences'])
        mutation_rate = float(mutation_rate)
        rec_rate = float(rec_rate)

        if ibd_len > 0:
            tmrca = (1 + k) / (ibd_len * 2 * (mutation_rate + rec_rate))
        else:
            tmrca = np.inf

        # Only consider pairs with recent shared ancestry
        if tmrca < tmrca_threshold:
            mut1 = mutation_status.get(i, 0)
            mut2 = mutation_status.get(j, 0)

            # How "recent" this ancestor is — closer to 1 if recent
            recency_weight = np.exp(-tmrca / scaling_factor)

            if mut1 == 1 and mut2 == 1:
                # Apply bonus: reduce distance
                tmrca *= (1 - recency_weight * (1 - bonus))
            elif mut1 != mut2:
                # Apply penalty: increase distance
                tmrca *= (1 + recency_weight * (penalty - 1))
            # else: both 0 → no change
        
        if not np.isfinite(tmrca):
            tmrca = 1e4

        distance_matrix[i, j] = tmrca
        distance_matrix[j, i] = tmrca  # symmetric
    
    # Condense and cluster
    condensed = squareform(distance_matrix)
    Z = linkage(condensed, method='average')
 
    return Z  

if __name__ == '__main__':
    # constant for all runs 
    mutation_position = 35001
    population_size = 1000
    sample_size = 100
    meta_path = "metas/phasedibd_meta.csv"
    vcf_path = sys.argv[1]

    # get all paths
    trees_path, sample_ID_path, zarr_path, dendogram_path, Z_output, mutation_rate, csv_output, csv_unique, rec_rate, index = get_paths(vcf_path)

    # use convert function
    cols, sweep, ht, pos = convert(zarr_path, sample_size, mutation_position)

    # run phaseibd analysis
    Z = run_phasedIBD(zarr_path, vcf_path, csv_output, csv_unique, ht, float(mutation_rate), float(rec_rate), sweep, pos, index)
    np.save(Z_output, Z)

    # run tmrca analysis to score and make dendogram
    run_tmrca_analysis(Z, cols, trees_path, sample_ID_path, dendogram_path, int(population_size), float(mutation_rate), meta_path, get_vcf_basename(vcf_path))

