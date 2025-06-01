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
from tmrca_phaseibd import get_paths, convert
from Z_analysis import get_Z_paths

MEOWSS_Z_paths = ["3/TMRCA/0.5/100/3_13_1000_0.00025_2.5e-05_0.5_100_0.96_20_100_Z.npy",
"3/TMRCA/0.7/100/3_13_1000_0.00025_2.5e-05_0.7_100_0.96_30_100_Z.npy",
"3/TMRCA/0.5/100/3_15_1000_0.00025_0.0025_0.5_100_0.96_30_100_Z.npy",
"3/TMRCA/0.7/100/3_15_1000_0.00025_0.0025_0.7_100_0.97_30_100_Z.npy",
"3/TMRCA/0.7/100/3_16_1000_0.0025_2.5e-05_0.7_100_0.97_30_100_Z.npy",
"3/TMRCA/0.7/100/3_17_1000_0.0025_0.00025_0.7_100_0.97_30_100_Z.npy",
"7/TMRCA/0.7/100/7_13_1000_0.00025_2.5e-05_0.7_100_0.97_20_100_Z.npy",]

for MEOWSS_Z_path in MEOWSS_Z_paths:
    parts = MEOWSS_Z_path.split("/")
    seed = parts[0]
    AF = parts[2]
    sample_size = parts[3]
    filename = parts[4]

    vcf_path = os.path.join(
        seed,
        "VCF_NEW",
        AF,
        sample_size,
        f"{filename[:-18]}"

    )


    meta_path = "metas/MEOWSS_meta.csv"
    sample_size = 100
    mutation_position = 35001
    population_size = 1000

    # get paths
    trees_path, sample_ID_path, zarr_path, dendogram_path, Z_output, mutation_rate, csv_output, csv_unique, rec_rate, index = get_paths(vcf_path) # uses with .vcf as input


    Z = np.load(MEOWSS_Z_path)

    # get cols
    cols, _, _, _ = convert(zarr_path, sample_size, mutation_position)

    # run tmrca analysis
    run_tmrca_analysis(
        Z, cols, trees_path, sample_ID_path, dendogram_path,
        int(population_size), float(mutation_rate), meta_path,
        get_vcf_basename(vcf_path) # uses with .vcf as input
        )

