from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist
import numpy as np
from scipy.stats import pearsonr
import pandas as pd
import os
import sys
import csv


def get_CC(Z1_path, Z2_path):
    """
    function that takes two Z values and gives the Cophenetic correlation between 2 dendograms (Z linkage matrix)

    """
    try:
        Z1 = np.load(Z1_path)
        Z2 = np.load(Z2_path)
        coph_dists1 = cophenet(Z1)
        coph_dists2 = cophenet(Z2)
        corr, _ = pearsonr(coph_dists1, coph_dists2)

        return corr

    except Exception as e:
        print(f"Error loading or comparing: {Z1_path} vs {Z2_path} -- {e}")
        
        return None


def get_Z_paths(vcf_path):
    """
    format of Z files
    get_vcfs_for_SINGER.txt --> 1/VCF_NEW/0.5/100/1_13_1000_0.00025_2.5e-05_0.5_100
    perfect: 1/PERFECT_TREES/0.5/Z/1_13_1000_0.00025_2.5e-05_0.5_100.npy
    singer: 1/SINGER/0.5/Z/1_13_1000_0.00025_2.5e-05_0.5_100_{index}.npy
    phasedibd: 1/PHASEIBD/0.5/Z/1_13_1000_0.00025_2.5e-05_0.5_100.npy
    MEOWSS : 1/TMRCA/0.5/100/1_17_1000_0.0025_0.00025_0.5_100_0.97_10_100_Z.npy

    """
    parts = vcf_path.split("/")
    seed = parts[0]
    AF = parts[2]
    sample_size = parts[3]
    filename = parts[4]
    params = filename.split("_")
    index = params[1]

    perfect_Z_path = os.path.join(
        seed,
        "PERFECT_TREES",
        AF,
        "Z",
        f"{filename}.npy"
    )

    # there are 100 singer Z files for each simulation!
    singer_Z_paths = []
    for i in range(0,100):
        singer_Z_path = os.path.join(
            seed,
            "SINGER",
            AF,
            "Z",
            f"{filename}_{i}.npy"
        )
        singer_Z_paths.append(singer_Z_path)


    phasedibd_Z_path = os.path.join(
        seed,
        "PHASEIBD",
        AF,
        "Z",
        f"{filename}.npy"
    )

    MEOWSS_Z_path = os.path.join(
        seed,
        "TMRCA",
        AF,
        sample_size,
        f"{filename}_0.97_10_100_Z.npy"
    )


    return perfect_Z_path, singer_Z_paths, phasedibd_Z_path, MEOWSS_Z_path, seed, index, AF

def average_out_singer(perfect_Z_path, singer_Z_paths):
    """
    Process SINGER Z files. Skip missing files, average available ones.
    """
    singer_CCs = []

    for path in singer_Z_paths:
        if os.path.exists(path):
            try:
                singer_CC = get_CC(perfect_Z_path, path)
                singer_CCs.append(singer_CC)
            except Exception as e:
                print(f"Error comparing {path}: {e}")
        else:
            print(f"Missing SINGER Z file: {path}")

    if len(singer_CCs) == 0:
        return None  # No data available

    return np.mean(singer_CCs)


def fill_row(perfect_Z_path, avg_singer_CC, phasedibd_Z_path, MEOWSS_Z_path, seed, index, AF):
    """
    get CC values for the rest and add them to a temp csv file which has to be merged later 

    """

    # calculate CCs
    phasedibd_CC = get_CC(perfect_Z_path, phasedibd_Z_path)
    MEOWSS_CC = get_CC(perfect_Z_path, MEOWSS_Z_path)

    row = [seed, index, AF, avg_singer_CC, phasedibd_CC, MEOWSS_CC] 
    csv_file = 'Z_analysis.csv'

    with open(csv_file, 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(row)
        
    return

if __name__ == '__main__':
    # load vcf file name
    vcf_path = sys.argv[1]

    # get info on paths
    perfect_Z_path, singer_Z_paths, phasedibd_Z_path, MEOWSS_Z_path, seed, index, AF = get_Z_paths(vcf_path)

    # process singer seperately as it has 100 values each simulation 
    avg_singer_CC = average_out_singer(perfect_Z_path, singer_Z_paths)

    # process the rest and fill rows 
    fill_row(perfect_Z_path, avg_singer_CC, phasedibd_Z_path, MEOWSS_Z_path, seed, index, AF)



