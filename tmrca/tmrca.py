# Import modules
import numpy as np
import zarr
import allel
import scipy.cluster.hierarchy as sch
import matplotlib
import matplotlib.pyplot as plt
import scipy.signal
from scipy.ndimage import gaussian_filter1d
import dask
from dask.base import compute
from itertools import combinations
import seaborn as sns
# import tskit
import pandas as pd
# meowcha added this below
from scipy.cluster.hierarchy import fcluster
from collections import Counter, defaultdict
import csv
import sys, os
import traceback


#prev student: Like Hamming distance code, this was also taken from Anushka Thawani. Adaptations were made to this on the 
#high-performance computer using shell script, but this could not be represented.

def convert(file, sample_size, mutation_position):
    '''
    This function extracts haplotypes sequences from a .vcz file 
    Adapted from: http://alimanfoo.github.io/2018/04/09/selecting-variants.html 
    
    Arguments:
        file: name of vcz file (converted from vcf from SLiM soft sweep simulation)
        
    Returns:
        ht: haplotype sequences for 200 individuals
        samp_freq: frequency of sweep mutation in sample
        cols: used to color dendrogram
        pos: position of variants 

    '''
    print("path is", file)

    # Open as a Dataset
    data = zarr.open_group(file, mode='r')
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
    
    return ht, pos, samp_freq, cols, sweep

def sliding_homozygosity(ht, pos, gts, genome_length, window):
    '''
    This function calculates the sliding homozygosity for all pairs of haplotypes.

    Arguments:
    gts : number of haplotype pairs
    ht : vector of haplotype sequences from convert() function in previous codeblock
    pos: position of variants from convert() function in previous codeblock

    Returns:
    homozygosities:  homozygosity of all haplotypes in ht in a array(?)
    '''
    # Initialise empty vector 
    hom = np.empty(shape=(genome_length,gts),dtype=np.float32)  

    reg =  slice(-100, -100, None)
    
    
    for x in range(0,genome_length):
        start  = x
        end = x + window
        try:    
            region = pos.locate_range(start,end)  
            # Check if the current window (region) is different from the previous window (reg) - makes code faster
            if region != reg:
                htx = ht[region]
                d = allel.pairwise_distance(htx, metric="hamming")
                hom[x,:] = 1-d  # 1-hamming distance = homozygosity
                reg = region
    
            else:
                hom[x,:] = 1-d

        except KeyError:
            pass
    
    return hom

def sliding_homozygosity_single(ht, pos, hap1_idx, hap2_idx, genome_length, window):
    '''
    Calculates sliding homozygosity between a single haplotype pair.

    Arguments:
        ht: haplotype sequences (2D array)
        pos: variant positions (SortedIndex)
        hap1_idx: index of first haplotype
        hap2_idx: index of second haplotype
        genome_length: length of genome
        window: size of sliding window

    Returns:
        hom: 1D array of homozygosity values
    '''
    hom = np.empty(genome_length, dtype=np.float32)
    reg = slice(-100, -100, None)

    for x in range(genome_length):
        start = x
        end = x + window
        try:
            region = pos.locate_range(start, end)
            if region != reg:
                htx = ht[region]
                hap1 = htx[:, hap1_idx]
                hap2 = htx[:, hap2_idx]
                # Calculate Hamming distance manually
                mismatches = np.sum(hap1 != hap2)
                homozygosity = 1 - mismatches / len(hap1) if len(hap1) > 0 else np.nan
                hom[x] = homozygosity
                reg = region
            else:
                hom[x] = homozygosity  # reuse last value if same window

        except KeyError:
            hom[x] = np.nan  # Or 0, depending on what makes more sense

    return hom


def sliding_homozygosity_single_slice(bp_start, bp_end, hom, ht, pos, hap1_idx, hap2_idx, window):
    '''
    Calculates sliding homozygosity between a single haplotype pair.

    Arguments:
        ht: haplotype sequences (2D array)
        pos: variant positions (SortedIndex)
        hap1_idx: index of first haplotype
        hap2_idx: index of second haplotype
        genome_length: length of genome
        window: size of sliding window

    Returns:
        hom: 1D array of homozygosity values
    '''
    reg = slice(-100, -100, None)
    print(f'calculating homozygosity for [{bp_start}, {bp_end})')

    for x in range(bp_start, bp_end):
        start = x
        end = x + window
        try:
            region = pos.locate_range(start, end)
            if region != reg:
                htx = ht[region]
                hap1 = htx[:, hap1_idx]
                hap2 = htx[:, hap2_idx]
                # Calculate Hamming distance manually
                mismatches = np.sum(hap1 != hap2)
                homozygosity = 1 - mismatches / len(hap1) if len(hap1) > 0 else np.nan
                hom[x] = homozygosity
                reg = region
            else:
                hom[x] = homozygosity  # reuse last value if same window

        except KeyError:
            hom[x] = np.nan  # Or 0, depending on what makes more sense

    return hom


def finding_troughs(smooth, mutation_position, threshold):
    '''
    This function finds troughs for a pair of haplotype sequences.

    Arguments:
        smooth: smoothed sliding window homozygosity for all pairs of sequences
        mutation_position: position of sweep mutation in genome (same as previous function)

    Returns:
        lower: position of breakpoint left of the sweep site
        upper: position of breakpoint right of the sweep site
        SHL: shared haplotype length
    '''

    #finding troughs
    try:
        troughs, _ = scipy.signal.find_peaks(-smooth) # - sign inverts the graph so the 'peaks' are our troughs
        troughs = troughs[smooth[troughs] < threshold]   # Extract troughs where homozygosity<threshold
    except IndexError as e:
        print("Index error in finding_troughs when extracting troughs", mutation_position)
        raise e

    #finds the peaks
    try:
        peaks, _ = scipy.signal.find_peaks(smooth)
    except IndexError as e:
        print("Index error in finding_troughs when finding peaks", mutation_position)
        raise e

    # Find positions of troughs flanking sweep site
    bp = np.searchsorted(troughs,mutation_position)  #search sorted finds index of position where mutation should be inserted in order to maintain the same order
    try:
        lower_trough = troughs[bp - 1] #index of sweep site -1
    except IndexError:
        print("Cannot find a trough on the left lower than threshold", threshold, "Setting lower_trough to 0")
        lower_trough = 0
    try:
        upper_trough = troughs[bp]  #index of sweep site
    except IndexError:
        print("Cannot find a trough on the right lower than threshold", threshold, "Setting upper_trough to", smooth.size - 1)
        upper_trough = smooth.size - 1

    # Find the average peak position around the sweep site
    try:
        highest = peaks[(peaks >= lower_trough) & (peaks <= upper_trough)]
        if highest.size != 0:
            highest = np.mean(highest)
        else:
            highest = (lower_trough+upper_trough)/2
    except IndexError as e:
        print("Index error in finding_troughs when finding avg peak position")
        raise e

    highest_float = highest
    highest = int(highest)

    # Get homozygosity for the highest peak, get whats half of it too
    # with np.printoptions(threshold=np.inf):
    #     print("smooth", smooth[0:10])
    highest_homozygosity = smooth[highest]
    half_homozygosity = highest_homozygosity / 2

    # LEFT SIDE — getting boundaries on lower_trough side
    try:
        left_seg = smooth[lower_trough:highest+1]
        if left_seg.size == 0: # ensure we don't attempt to find argmin of an empty sequence
            return 0, 0, 0
        left_idx = np.where((left_seg[:-1] < half_homozygosity) & (left_seg[1:] >= half_homozygosity))[0]
        if left_idx.size > 0:
            lower = lower_trough + left_idx[0]
        else:
            lower = lower_trough + np.argmin(np.abs(left_seg - half_homozygosity))
    except IndexError as e:
        print("IndexError in LEFT SHL boundary", traceback.format_exc())
        raise e

    # RIGHT SIDE — getting boundaries on upper_trough side
    try:
        right_seg = smooth[highest:upper_trough+1]
        if right_seg.size == 0:
            return 0, 0, 0
        right_idx = np.where((right_seg[:-1] >= half_homozygosity) & (right_seg[1:] < half_homozygosity))[0]
        if right_idx.size > 0:
            upper = highest + right_idx[0] + 1
        else:
            upper = highest + np.argmin(np.abs(right_seg - half_homozygosity))
    except IndexError as e:
        print("IndexError in RIGHT SHL boundary", traceback.format_exc())
        raise e

    SHL = upper - lower
    print("Lower", lower, "Upper", upper, "SHL", SHL,
        "Highest", highest_float,
        "Old lower", (lower_trough+highest_float)/2,
        "Old upper", (upper_trough+highest_float)/2,
        "Old SHL", (upper_trough+highest_float)/2 - (lower_trough+highest_float)/2 )

    return int(lower), int(upper), SHL


## function using each haplotype pair to return SHL and find lower and upper limits of SHL
# uses finding_troughs()
def find_breakpoint(haplotype_pair, mutation_position, points, threshold):
    '''
    For a pair of sequences, this function smoothes the sliding homozygosity and returns the SHL
    Arguments:haplotype_pair
        haplotype_pair: a pair of haplotype sequences
        
    Returns:
        lower: position of breakpoint left of the sweep site
        upper: position of breakpoint right of the sweep site
        SHL: shared haplotype length
    '''
    
    smooth = gaussian_filter1d(haplotype_pair, points)
    try:
        lower, upper, SHL = finding_troughs(smooth, mutation_position, threshold)
    except IndexError:
        print("Index error in find_breakpoint!!", traceback.format_exc())
        lower = -1.3
        upper = -1.3
        SHL = -1.3
        
    return lower, upper, SHL, smooth

def find_lower_and_upper_troughs_fast(ht, pos, hap1_idx, hap2_idx, genome_length, mutation_position, points, threshold, window, check_interval):
    '''
    For a pair of sequences, this function smoothes the sliding homozygosity and returns the SHL
    Arguments:
        check_interval: The interval at which we sample homozygosity before checking for trough.
        
    Returns:
        lower: position of breakpoint left of the sweep site
        upper: position of breakpoint right of the sweep site
        SHL: shared haplotype length
    '''
    blur_sd_threshold = 4
    blur_radius = points * blur_sd_threshold
    blur_window_size = 2 * blur_radius + 1
    middle = mutation_position
    gaussian_kernel = scipy.signal.windows.gaussian(blur_window_size, std=points)
    gaussian_kernel /= gaussian_kernel.sum() # normalize kernel

    # Set default homozygosity to a crazy value to make it obvious when we access a yet undefined homozygosity
    homozygosities = np.full(genome_length, -99999.0)
    smooth_homozygosities = np.full(genome_length, -99999.0)

    # The lowest-index bp for which we calculated homozygosity.
    min_discovered = middle - blur_radius
    # The highest-index bp for which we calculated homozygosity.
    max_discovered = middle + blur_radius
    
    # Calculate homozygosity for the initial slice.
    # Initially, we calculate for the area (middle - blur_radius, middle + blur_radius)
    # so that values around the middle smooth correctly.
    sliding_homozygosity_single_slice(min_discovered - 1, max_discovered + 1,
                                      homozygosities, ht, pos, hap1_idx,
                                      hap2_idx, window)

    lower_trough = 0
    upper_trough = genome_length
    
    # Find closest left trough by checking homozygosity for check_interval bps at a time.
    while min_discovered >= 0:
        next_min_discovered = max(0, min_discovered - check_interval)
        sliding_homozygosity_single_slice(next_min_discovered, min_discovered,
                                      homozygosities, ht, pos, hap1_idx,
                                      hap2_idx, window)
        smooth_part = scipy.signal.convolve(homozygosities[next_min_discovered:min_discovered+blur_window_size],
                                            gaussian_kernel, mode='valid')

        min_smooth = next_min_discovered + blur_radius # index of the smallest value that has been smoothed out
        smooth_homozygosities[min_smooth:min_smooth+smooth_part.size] = smooth_part
    
        #print(smooth_part)
        troughs, _ = scipy.signal.find_peaks(-smooth_part) # invert sign to find troughs
        troughs = troughs[smooth_part[troughs] < threshold] # only allow troughs below threshold
        troughs = [min_smooth + trough for trough in troughs] # map values in troughs to BP

        min_discovered = next_min_discovered

        if troughs: # If this segment contains a trough
            lower_trough = troughs[-1] # take the highest index trough (scipy's find_peaks return sorted ascending array) 
            break

    # Find closest right trough by checking homozygosity for check_interval bps at a time.
    while max_discovered >= 0:
        next_max_discovered = min(genome_length, max_discovered + check_interval)
        sliding_homozygosity_single_slice(max_discovered, next_max_discovered,
                                      homozygosities, ht, pos, hap1_idx,
                                      hap2_idx, window)
        smooth_part = scipy.signal.convolve(homozygosities[max_discovered-blur_window_size:next_max_discovered],
                                            gaussian_kernel, mode='valid')
    
        min_smooth = max_discovered-blur_window_size+blur_radius # index of the smallest value that has been smoothed out
        smooth_homozygosities[min_smooth:min_smooth+smooth_part.size] = smooth_part
        
        troughs, _ = scipy.signal.find_peaks(-smooth_part) # invert sign to find troughs
        troughs = troughs[smooth_part[troughs] < threshold] # only allow troughs below threshold
        troughs = [min_smooth + trough for trough in troughs] # map values in troughs to BP

        max_discovered = next_max_discovered

        if troughs: # If this segment contains a trough
            upper_trough = troughs[0] # take the lowest index trough (scipy's find_peaks return sorted ascending array) 
            break

    print(lower_trough, upper_trough)

    return lower_trough, upper_trough, upper_trough - lower_trough, smooth_homozygosities

# Calculating Kij (No. of SNPs in each SHL)
def find_snp(n,gts,ht,results_computed_1,pos):
    '''
    This function finds the number of SNPs over the shared haplotype length for all pairs of haplotype sequences
    
    Arguments:
        sample_size: sample size to calculate number of haplotype sequences
        gts: number of haplotype pairs
        ht: haplotype sequnces
        results_computed_1: output from find_breakpoint function
        pos: position of variants from convert() function 
        
    Returns:
        diffs: number of SNP differences for all pairs of haplotype sequences
        
    '''
    n = sample_size * 2
    pairwise = []
    for combo in combinations(list(range(0,n)), 2): 
        pairwise.append(combo)

    diffs = np.empty(shape=(gts),dtype=np.float32)
    for i in range(gts):
        pair = ht[:,pairwise[i]]
        try:
            start = results_computed_1[i,1]
            stop = results_computed_1[i,2]

            window_pos = pos.locate_range(start, stop)
            window = pair[window_pos]

            d = allel.pairwise_distance(window, metric = "hamming")

            diffs[i]=d 

        except KeyError:
            diffs[i]=0 # set it to 0 bcs having a lot of negative values broke stuff
            print("key error in find_snp", traceback.format_exc())
    return diffs


#Calculating Tij, Time to Common Ancestor (TMRCA) using mutation rate and number of SNPs
# 𝜏_𝑖𝑗=(𝑘_𝑖𝑗+1)/(2ℓ_𝑖𝑗 (𝑟+𝜇))
#for array number x:
    #read in vcf file
    # read in population size, mu, r from population parameters csv file


def calculating_Tij(sample_size, mutation_position, genome_length, window, points, threshold):
    '''
    This function calculates Tij, Time to Common Ancestor (TMRCA) using mutation rate and number of SNPs for each haplotype pair.
    It uses the convert() to convert files from vcf 
    It uses the sliding_heterozygosity() to calculate heterozygosity in the window for all pairs of haplotypes

    
    Arguments:
        array_index: array index or combination number of simulation
        seed: seed number of simulation
        
        genome_length: length of genome (in SLiM simulation)          
        ht : haplotypes (what variable structure?)

        window: length of sliding window
        threshold: threshold above which troughs are ignored
        points: number of points to use for 1D-gaussian filter (see scipy documentation)
        
    Returns:
        Tij: Time to most recent Common Ancestor (TMRCA)
        cols: used to color dendrogram
  
    '''    

    # read in population size, mu, r from population parameters csv

    # Extract haplotype sequences from .vcf file
    ht, pos, samp_freq, cols, sweep = convert(filename, sample_size, mutation_position)

    # Calculate sliding homozygosity for all pairs of haplotype sequences
    no_haplotypes = sample_size * 2
    gts = int((no_haplotypes*(no_haplotypes-1))/2)
    homozygosities = sliding_homozygosity(ht, pos, gts, genome_length, window)

    print("homozygosities:", homozygosities, "len", len(homozygosities))

    
    # Find SHL for all pairs of haplotype sequences 
    hom_dask = dask.array.from_array(homozygosities, chunks=(genome_length, 1)) # type: ignore # creates a dask array
    homozygosities = []
    results = dask.array.apply_along_axis(find_breakpoint, 0, hom_dask, mutation_position, points, threshold) # type: ignore #applies find_breakpoint() along the array
    results_computed = results.compute()

    # Manipulating the dataframe to make it easier to process
    results_computed = np.transpose(results_computed)
    index = np.asarray(range(0,gts))
    index = np.expand_dims(index, axis=0)
    results_computed_1 = np.concatenate((index.T, results_computed), axis=1)
    
    
    # Calculate the TMRCA from the SHLs and number of SNPs
    shls = results_computed_1[:,3]   # SHLs for all pairs of haplotype sequences 
    shls[shls<=0] = genome_length
    diffs = find_snp(sample_size, gts, ht, results_computed_1, pos)
    Tij = (1+diffs)/(2*shls*(recombination_rate + mutation_rate)) # TMRCA metric for all pairs of haplotype sequences 

    
    # Remove negative and non-integer TMRCA values
    impute = np.nanmean(Tij)        #impute is the mean of SNP array without any NAN values
    print("impute", impute)
    x = np.isfinite(Tij)            #x is a boolean mask array of only finite values (cannot be infinite or NAN)
    for i in np.where(x == 0)[0]:   #for all indices where there is a non-finite number:
        Tij[i] = impute             #replace NaN with inpute value
    Tij[Tij<=0] = impute            # replace all negative numbers wih inpute

    return Tij, cols


def analysis(sample_size, mutation_position, genome_length, window, points, threshold): 

    '''
    This function plots a dendrogram and colours red the tips that have the sweep mutations. 
    It uses calculating_Tij().
    Only the output vcf from the SLIM simulation is input.
    global variables genome_length, window, threshold, points being used

    Arguments:
    
    global variables/ variables from sub functions:
    window: length of sliding window
    threshold: threshold above which troughs are ignored
    points: number of points to use for 1D-gaussian filter (see scipy documentation)
    cols: colours haplotype branches red if they have sweep mutation. From convert().



    Returns:
    output dendrogram in pdf
    '''
    #all convert(), etc etc to end up with tij
    Tij, cols = calculating_Tij(sample_size, mutation_position, genome_length, window, points, threshold)
    
    # Hierachical Clustering, store in Z
    Z = sch.linkage(Tij, method = 'average') #why do we use the Farthest Point Algorithm? changed to average (UPGMA) , check after

    ## Plot dendrogram without colouring branches
    # updating matplotlib font settings
    matplotlib.rcParams.update({'font.size': 24})
    fig = plt.figure(figsize=(30, 12))
    gs = matplotlib.gridspec.GridSpec(2, 1, hspace=0.1, wspace=1, height_ratios=(1,1)) # type: ignore

    ax_dend = fig.add_subplot(gs[0, 0])
    sns.despine(ax=ax_dend, offset=5, bottom=True, top=True)
    dd = sch.dendrogram(Z,color_threshold=0,above_threshold_color='#808080',ax=ax_dend) # if above colour threshold, set colour to grey 

    ls = []
    for leaf, leaf_color in zip(plt.gca().get_xticklabels(), dd["leaves_color_list"]):   #leaves_color_list is A list of color names. The k’th element represents the color of the k’th leaf.
        leaf.set_color(cols[int(leaf.get_text())])
        ls.append(int(leaf.get_text()))

    ax_dend.set_ylabel('Haplotype age/generations',fontsize=24)
    ax_dend.set_title('Haplotype clusters',fontsize=24)



    # Plot dendrogram and colour branches
    ax_dend_2 = fig.add_subplot(gs[1, 0])
    
    dflt_col = "#808080"
    
    link_cols = {}
    for i, i12 in enumerate(Z[:,:2].astype(int)):
        c1, c2 = (link_cols[x] if x > len(Z) else cols[x] for x in i12)
        link_cols[i+1+len(Z)] = c1 if c1 == c2 else dflt_col

    sns.despine(ax=ax_dend_2, offset=5, bottom=True, top=True)
    dd = sch.dendrogram(Z,link_color_func=lambda x: link_cols[x], ax=ax_dend_2)

    ls = []
    for leaf, leaf_color in zip(plt.gca().get_xticklabels(), dd["leaves_color_list"]):
        leaf.set_color(cols[int(leaf.get_text())])
        ls.append(int(leaf.get_text()))

    ax_dend_2.set_ylabel('Haplotype age/generations',fontsize=24)
    
    
    # Save dendrogram
    output_directory = f'{seed_index}/TMRCA/{frequency}/{sample_size}' 
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    pdf_output = f'{output_directory}/{basename}_update1.pdf'
    Z_output = f'{output_directory}/{basename}_update1_Z.npy'
    cols_output = f'{output_directory}/{basename}_update1_cols.npy'
    print (pdf_output)
    print(Z_output)
    print(cols_output)
    plt.savefig(pdf_output)
    np.save(Z_output, Z)
    np.save(cols_output, cols)
    

if __name__ == '__main__':
        
    # These are the same for everything
    genome_length = 70001
    mutation_position = int((genome_length+1)/2)


    threshold = 0.87  # setting threshold to calculate trough points in each haplotype
    points=280 # smoothing 
    window=600 # sliding window for homozygosity 

    # These are different for every vcf file
    try:
        filename = sys.argv[1]
        threshold = float(sys.argv[2])
        points = int(sys.argv[3])
        window = int(sys.argv[4])
    except IndexError:
        print("Usage: python3 tmrca.py <vcz_filename> <threshold> <points> <window>")
        exit()

    basename = os.path.basename(filename).strip()[:-4] # removes the .vcz extension from filename
    basename_components = basename.split('_')
    basename += f'_{threshold}_{points}_{window}'
    global seed_index
    seed_index = int(basename_components[0])
    global population_size
    population_size = int(basename_components[2])
    global mutation_rate
    mutation_rate = float(basename_components[3])
    global recombination_rate
    recombination_rate = float(basename_components[4])
    global frequency
    frequency = float(basename_components[5])
    sample_size = int(basename_components[6])

    print("population_size", population_size)
    print("mutation_rate", mutation_rate)
    print("recombination_rate", recombination_rate)
    print("sample_size", sample_size)
    analysis(sample_size, mutation_position, genome_length, window, points, threshold)