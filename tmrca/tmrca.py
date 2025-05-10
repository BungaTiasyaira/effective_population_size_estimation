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
import dask.bag as db
from dask.base import compute
from itertools import combinations
from functools import partial
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
    # print(f'calculating homozygosity for [{bp_start}, {bp_end}]')

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



def find_lower_and_upper_troughs_fast(ht, pos, hap1_idx, hap2_idx, genome_length, mutation_position, points, threshold, window, check_interval):
    '''
    For a pair of sequences, this function smoothes the sliding homozygosity and returns the SHL
    The sliding window begins at the 35000 (location of mutation), moves left and right, looks ahead "check_interval" base pairs
    calculates homozygosity using sliding_homozygosity_single_slice function and does gaussian smoothing.
    Checks for troughs below the set threshold, stops when found on both left and right, gets SHL. 
    
    Arguments:
		ht: haplotype sequences (2D array)
        pos: variant positions (SortedIndex)
        hap1_idx: index of first haplotype
        hap2_idx: index of second haplotype
        genome_length: length of genome set globally 
        mutation_position: 35000 or middle of genome
        points: smoothing points set globally 
        threshold: trough threshhold set globally 
        window: size of sliding window
        check_interval: The interval at which we sample homozygosity before checking for trough.
        
    Returns:
        lower: position of breakpoint left of the sweep site
        upper: position of breakpoint right of the sweep site
        smooth_homozygosities: array of smoothed homozygosities
    '''
    # print idx 
    print(hap1_idx, hap2_idx)

    # Set parameters for Gaussian smoothing
    blur_sd_threshold = 4
    blur_radius = points * blur_sd_threshold
    blur_window_size = 2 * blur_radius + 1
    middle = mutation_position
    gaussian_kernel = scipy.signal.windows.gaussian(blur_window_size, std=points)
    gaussian_kernel /= gaussian_kernel.sum() # normalize kernel

    # Normalized Gaussian kernel used for smoothing homozygosity values
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

	# Defaults troughs to entire genome length
	# will be updated when real troughs are found.
    lower_trough = 0
    upper_trough = genome_length - 1
    lower_peak = 0
    upper_peak = genome_length - 1

    # Defaults peaks to middle, will update when peaks found 
    lower_peak = middle - 1
    upper_peak = middle + 1


    def handle_smooth_part(homozygosities_to_scan, min_smooth, direction):
        '''
	    Smooths a region of the homozygosity array with a Gaussian.
	    Finds peaks (maxima) above the threshold.
	    Returns positions adjusted to genome coordinates.
	    
	    Arguments:
		    homozygosities_to_scan: 1D window of homozygosity values.
		    min_smooth:the starting index in the genome to smooth 
            direction: left or right
		    
		  Returns:
		    final_peak: peak identified  
            final_trough: trough identified 
	    '''
        # Applies gaussian smoothing to homozygosities_to_scan
        smooth_part = scipy.signal.convolve(homozygosities_to_scan,
                                            gaussian_kernel, mode='valid') # mode='valid': results in a shorter output than the input. 

		# stores the smooth_part back into global smooth_homozygosities array, aligned to genomic coordinates
        smooth_homozygosities[min_smooth:min_smooth+smooth_part.size] = smooth_part
    
        # find peaks and troughs
        peaks, _ = scipy.signal.find_peaks(smooth_part) # finding peak
        troughs, _ = scipy.signal.find_peaks(-smooth_part) # finding trough, inverse peak

        # Adjusts peak indices to genome-wide coordinates by adding the min_smooth offset.
        peaks = [min_smooth + peak for peak in peaks] # map values in peaks to BP
        troughs = [min_smooth + trough for trough in troughs]

        def get_peak_trough(peaks, troughs, direction):
            '''
            Based on the identified peaks and troughs, it finds a pair of peak and trough iteratively which has a sufficient
            difference in homozygosity values. The pair will be returned as the final peak and trough 

            Arguments:
                peaks: list of peaks identified
                troughs: list of troughs identified 
                direction: left or right 

            Returns:
                final_peak: peak identified  
                final_trough: trough identified 

            '''
            mutation_pos = 35000   
            # Handle empty peak/trough cases
            if not peaks or not troughs:
                # print('peak or trough list is empty')
                return None, None
    
            # Initialize direction
            if direction == 'left':
                peaks = sorted(peaks, reverse=True)
                troughs = sorted(troughs, reverse=True)
            elif direction == 'right':
                peaks = sorted(peaks)
                troughs = sorted(troughs)

            # Walk through peaks and try to find a matching trough
            for i, peak in enumerate(peaks):
                H_peak = smooth_homozygosities[peak]
                trough_threshold = H_peak * threshold

                # Find the first trough that comes **after** the peak in the direction
                for trough in troughs:
                    if (direction == 'left' and trough < peak) or (direction == 'right' and trough > peak):
                        H_trough = smooth_homozygosities[trough]
                        if H_trough < trough_threshold:
                            # print(direction, "peak", peak, "trough", trough)
                            return peak, trough  # Found a valid pair
                        else:
                            break  # Trough is too high — try next peak
            print("no valid peak-trough found")
            return None, None  # No valid peak-trough pair found

        final_peak, final_trough = get_peak_trough(peaks, troughs, direction)
 
        # Returns the detected peak and trough positions.
        return final_peak, final_trough 

    # Find closest left trough by checking homozygosity for check_interval bps at a time.
    while min_discovered > 0:
        # Computes next region to evaluate on the left.
        next_min_discovered = max(0, min_discovered - check_interval)
        # Fills in homozygosity values for this region.
        sliding_homozygosity_single_slice(next_min_discovered, min_discovered,
                                      homozygosities, ht, pos, hap1_idx,
                                      hap2_idx, window)
        
    
        # Prepares a slice for smoothing and its minimum smoothed coordinate.
        homozygosities_to_scan = homozygosities[next_min_discovered:min_discovered+blur_window_size]        
        min_smooth = next_min_discovered + blur_radius # index of the smallest value that will be smoothed out

        # find final peak and trough after smoothing using 'handle_smooth_part' function with 'get_peak_trough' 
        final_peak, final_trough = handle_smooth_part(homozygosities_to_scan, min_smooth, direction = 'left')

        min_discovered = next_min_discovered

        # If final trough found, asign it to lower_trough and break 
        if final_trough is not None: 
            lower_trough = final_trough
            lower_peak = final_peak  
            break
    
    if lower_trough == 0: # If still not found
        homozygosities_to_scan = np.concatenate((
            homozygosities[0:blur_radius],
            homozygosities[0:blur_window_size]
        ))
        min_smooth = 0
        final_peak, final_trough = handle_smooth_part(homozygosities_to_scan, min_smooth, direction = 'left')
        if final_trough:
            lower_trough = final_trough
            lower_peak = final_peak


    # Find closest right trough by checking homozygosity for check_interval bps at a time.
    while max_discovered < genome_length:
        # Computes next region to evaluate on the right.
        next_max_discovered = min(genome_length, max_discovered + check_interval)
        
        # Fills in homozygosity values for this region.
        sliding_homozygosity_single_slice(max_discovered, next_max_discovered,
                                      homozygosities, ht, pos, hap1_idx,
                                      hap2_idx, window)
        
        # Prepares a slice for smoothing and its minimum smoothed coordinate.
        homozygosities_to_scan = homozygosities[max_discovered-blur_window_size:next_max_discovered]
        min_smooth = max_discovered-blur_window_size+blur_radius # index of the smallest value that has been smoothed out
        
        # find final peak and trough after smoothing using 'handle_smooth_part' function with 'get_peak_trough' 
        final_peak, final_trough = handle_smooth_part(homozygosities_to_scan, min_smooth, direction = 'right')

        max_discovered = next_max_discovered

        # If final trough found, asign it to upper_trough and break 
        if final_trough is not None: 
            upper_trough = final_trough 
            upper_peak = final_peak
            break

    if upper_trough == genome_length - 1: # If still not found
        homozygosities_to_scan = np.concatenate((
            homozygosities[genome_length - blur_window_size:genome_length],
            homozygosities[genome_length - 1:genome_length - blur_radius:-1]
        ))
        min_smooth = genome_length - blur_radius
        final_peak, final_trough = handle_smooth_part(homozygosities_to_scan, min_smooth, direction = 'right')
        if final_trough:
            upper_trough = final_trough
            upper_peak = final_peak

    # get all troughs available in smooth_homozygosities
    all_troughs, _ = scipy.signal.find_peaks(-smooth_homozygosities)

    # calculating distance to mutation 
    distance_to_mut_lower_peak = abs(mutation_position - lower_peak)
    distance_to_mut_upper_peak = abs(mutation_position - upper_peak)
    difference_in_peaks = abs(smooth_homozygosities[lower_peak] - smooth_homozygosities[upper_peak])

    # if difference in peak is insignificant, I will pick the one closer to the mutation
    # if difference in peak is significant, I will rely on height, however if theyre far away, I will go with the closer
    def choose():
        def choose_closer():
            print('Therefore, choosing the closer one.')
            return "left" if distance_to_mut_lower_peak < distance_to_mut_upper_peak else "right" 
        def choose_higher():
            print('Therefore, choosing the higher one.')
            return "left" if smooth_homozygosities[lower_peak] > smooth_homozygosities[upper_peak] else "right"

        # print('Deciding which peak to pick...')
        print(f'Lower peak - {lower_peak}, {smooth_homozygosities[lower_peak]:.3}')
        print(f'Upper peak - {upper_peak}, {smooth_homozygosities[upper_peak]:.3}')
        print(f'Difference: {difference_in_peaks:.3}') # :.3 means to 3 significant figures


        # if difference in peak is insignificant, or they are far away, I will pick the one closer to the mutation
        if difference_in_peaks < 0.01:
            print(f'Difference in peaks is insignificant')
            return choose_closer()

        # if one of them is far away, pick the closer one 
        if distance_to_mut_lower_peak > 5000 or distance_to_mut_upper_peak > 5000:
            print('One of them is far away...')
            return choose_closer()
        
        # Otherwise, rely on height
        print(f'Both peaks are close, difference in peaks is significant...')
        return choose_higher()

    decision = choose()
    print(f'Chose {decision}.')
    assert(decision == "left" or decision == "right")
        
    if decision == "left":
        selected_threshold = smooth_homozygosities[lower_peak]*0.87
        chosen_peak = lower_peak
        notchosen_peak = upper_peak
        # scan right for valid troughs
        for trough in sorted(all_troughs):
            if trough > mutation_position and smooth_homozygosities[trough] < selected_threshold:
                upper_trough = trough
                break
    elif decision == "right":
        selected_threshold = smooth_homozygosities[upper_peak]*0.87
        chosen_peak = upper_peak
        notchosen_peak = lower_peak
        # scan left for valid troughs
        for trough in sorted(all_troughs, reverse=True):
            if trough < mutation_position and smooth_homozygosities[trough] < selected_threshold:
                lower_trough = trough
                break

    return lower_trough, upper_trough, chosen_peak, notchosen_peak, smooth_homozygosities

def calculate_SHL(lower, upper, smooth_homozygosities):
    peaks, _ = scipy.signal.find_peaks(smooth_homozygosities[lower:upper])
    highest = lower + round(np.mean(peaks) if peaks.size > 0 else (lower + upper) / 2)

    trough_average = (smooth_homozygosities[lower] + smooth_homozygosities[upper]) / 2
    half_threshold = (smooth_homozygosities[highest] - trough_average) / 2 + trough_average

    # Find first and last value greater than half_threshold
    # If no such value is found, left will be equal to lower
    # and right will be equal to upper
    left = np.argmax(smooth_homozygosities[lower:upper] >= half_threshold) + lower
    right = (upper - 1) - np.argmax(smooth_homozygosities[upper-1::-1] >= half_threshold)

    SHL = right - left

    return left, right, SHL

def process_haplotype_pair(pair_of_haplotypes, ht, pos, genome_length, mutation_position, points, threshold, window, check_interval):
    h1, h2 = pair_of_haplotypes
    lower_trough, upper_trough, _, _, smooth_homozygosities = find_lower_and_upper_troughs_fast(ht, pos, h1, h2, genome_length, mutation_position, points, threshold, window, check_interval)
        
    left, right, SHL = calculate_SHL(lower_trough, upper_trough, smooth_homozygosities)

    # now i want to immediately calculate SNP (Kij)
    try:
        region = pos.locate_range(left, right)
        htx = ht[region]
        hap1 = htx[:, h1]
        hap2 = htx[:, h2]
        # Calculate Hamming distance manually aka SNP diffs
        diff = np.sum(hap1 != hap2)        
    except KeyError:
        diff = 0 # set it to 0 bcs having a lot of negative values broke stuff
        print("key error in process_haplotype_pair (set to 0) for this below:")
        print("h1: " + h1, "h2: " + h2)
        print("left: " + left, "right: " + right, "SHL: " + SHL)
        print(traceback.format_exc())

    print(h1, h2)
    print("Left:", left, "Right:", right, "SHL:", SHL, "Diff:", diff)
    
    return [SHL, diff]

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

    # Extract haplotype sequences from .vcz file
    ht, pos, samp_freq, cols, sweep = convert(filename, sample_size, mutation_position)

    # get pairs (non repeating) of haplotypes
    # Calculate sliding homozygosity for all pairs of haplotype sequences using fast new method
    haplotype_combinations = list(combinations(range(ht.shape[1]), 2))
    process_haplotype_pair_func = partial(
        process_haplotype_pair,
        ht=ht,
        pos=pos,
        genome_length=genome_length,
        mutation_position=mutation_position,
        points=points,
        threshold=threshold,
        window=window,
        check_interval=check_interval
    )
    bag = db.from_sequence(haplotype_combinations, npartitions=4)
    results_array = bag.map(process_haplotype_pair_func).compute()
    
    print('finished calculating SHLs and diffs!')
    
    # turn into 1D numpy array
    results_array = np.array(results_array) # shape: (19900, 2)

    # Extract SHLs and diffs (SNPs)
    shls = results_array[:, 0]
    diffs = results_array[:, 1]

    # Fix SHLs if needed
    shls[shls <= 0] = genome_length

    # Calculate the TMRCA from the SHLs and number of SNPs
    Tij = (1+diffs)/(2*shls*(recombination_rate + mutation_rate)) # TMRCA metric for all pairs of haplotype sequences 
    print("calculated Tij!")
    
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
    Z = sch.linkage(Tij, method = 'average') # use average (UPGMA)

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

    pdf_output = f'{output_directory}/{basename}_updateMANTKU.pdf'
    Z_output = f'{output_directory}/{basename}_updateMANTKU_Z.npy'
    cols_output = f'{output_directory}/{basename}_updateMANTKU_cols.npy'
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


    # threshold = 0.87  # setting threshold to calculate trough points in each haplotype
    # points=280 # smoothing 
    # window=600 # sliding window for homozygosity 

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
    no_haplotypes = sample_size*2
    global gts
    gts = int((no_haplotypes*(no_haplotypes-1))/2)
    global check_interval
    check_interval = 1000


    print("population_size", population_size)
    print("mutation_rate", mutation_rate)
    print("recombination_rate", recombination_rate)
    print("sample_size", sample_size)

    analysis(sample_size, mutation_position, genome_length, window, points, threshold)