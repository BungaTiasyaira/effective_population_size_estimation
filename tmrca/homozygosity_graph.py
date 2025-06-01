import matplotlib.pyplot as plt
from tmrca import find_lower_and_upper_troughs_fast, convert, sliding_homozygosity_single
import numpy as np
import math

# setting parameters manually 
sample_size = 100
genome_length = 70001
window = 100
check_interval = 5000
mutation_position = int((genome_length+1)/2)
ht, pos, samp_freq, cols, sweep = convert('1_investigating_100mutants/ZARR/0.5/100/1_14_1000_0.00025_0.00025_0.5_100.vcz', sample_size, mutation_position)
no_haplotypes = sample_size * 2
gts = int((no_haplotypes*(no_haplotypes-1))/2)
points = 50
threshold = 0.87

# run function
lower_trough, upper_trough, chosen_peak, notchosen_peak, smooth_homozygosities = find_lower_and_upper_troughs_fast(ht, pos, 6, 63, genome_length, mutation_position, points, threshold, window, check_interval)

# print final values found 
print('lower  trough:', lower_trough, f"{smooth_homozygosities[lower_trough]:.2f}",
      'peak:', chosen_peak, f"{smooth_homozygosities[chosen_peak]:.2f}")

print('upper  trough:', upper_trough, f"{smooth_homozygosities[upper_trough]:.2f}",
      'peak:', chosen_peak, f"{smooth_homozygosities[chosen_peak]:.2f}")


# setting x values 
xmin = np.argmax(smooth_homozygosities>=0)
xmax = xmin + (smooth_homozygosities[smooth_homozygosities>=0]).shape[0]
selected_homozygosities = smooth_homozygosities[xmin:xmax]
x_values = list(range(xmin, xmin+len(selected_homozygosities)))

# plot smooth homozygosities
plt.plot(x_values, selected_homozygosities)

# plot points for troughs
plt.plot(lower_trough, smooth_homozygosities[lower_trough], marker='o', color='red', markersize=7)
plt.plot(upper_trough, smooth_homozygosities[upper_trough], marker='o', color='red', markersize=7)

# plot points for peaks 
plt.plot(chosen_peak, smooth_homozygosities[chosen_peak], marker='o', color='green', markersize=7)
plt.plot(notchosen_peak, smooth_homozygosities[notchosen_peak], marker='o', color='orange', markersize=7)

plt.plot(35000, smooth_homozygosities[35000], marker='*', color='purple', markersize=7)
plt.xlabel("Genomic Position")
plt.ylabel("Homozygosity")
plt.title("Homozygosity (fast)")
plt.savefig("tmrca/homozygosity.png")  # Saves the plot to a PNG file
plt.close()  # Closes the plot so it doesn't display in interactive environments
