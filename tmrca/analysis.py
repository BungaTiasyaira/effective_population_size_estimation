import numpy as np
import scipy.cluster.hierarchy as sch
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from analysis_trees import get_mapping_of_haplotype_to_generation
import sys, os

# Kit this file is just the analysis part from tmrca.py I basically copied it here
def load_data_for_analysis(Z_path, cols_path):
    Z = np.load(Z_path, allow_pickle=True)
    cols = np.load(cols_path, allow_pickle=True).item() # .item() is needed to load as dict
    return Z, cols

try:
    vcz_filename = sys.argv[1]
    threshold = float(sys.argv[2])
    points = int(sys.argv[3])
    window = int(sys.argv[4])
    basename = os.path.basename(vcz_filename)[:-4]
except IndexError:
    print("Usage: python3 analysis.py <vcz_filename> <threshold> <points> <window>")
    exit()

def get_paths(basename, threshold, points, window):
    basename_components = basename.split('_')
    seed_index = int(basename_components[0])
    burnin_number = int(basename_components[1])
    population_size = int(basename_components[2])
    mutation_rate = float(basename_components[3])
    recombination_rate = float(basename_components[4])
    frequency = float(basename_components[5])
    sample_size = int(basename_components[6])

    file_general_format = '_'.join(basename_components[:-1]) # exclude sample_size from filename

    tmrca_folder = f'{seed_index}/TMRCA/{frequency}/{sample_size}'
    tmrca_output_format = f'{tmrca_folder}/{file_general_format}_{sample_size}_{threshold}_{points}_{window}'
    Z_path = f'{tmrca_output_format}_update1_Z.npy'
    cols_path = f'{tmrca_output_format}_update1_cols.npy'

    trees_folder = f'{seed_index}/TREES'
    trees_path = f'{trees_folder}/{file_general_format}.trees'

    sample_ID_folder = f'{seed_index}/sample_ID/{frequency}/{sample_size}'
    sample_ID_path = f'{sample_ID_folder}/{file_general_format}_{sample_size}.txt'

    output_folder = f'{seed_index}/TMRCA_dendograms/{frequency}/{sample_size}'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    output_path = f'{output_folder}/{file_general_format}_{sample_size}_{threshold}_{points}_{window}_update1_analysis.pdf'

    return Z_path, cols_path, trees_path, sample_ID_path, output_path

Z_path, cols_path, trees_path, sample_ID_path, output_path = get_paths(basename, threshold, points, window)

Z, cols = load_data_for_analysis(Z_path, cols_path)
haplotype_to_generation = get_mapping_of_haplotype_to_generation(trees_path,
                                                                 sample_ID_path)

generation_colormap = matplotlib.colormaps['gist_rainbow']
generation_color = {}
list_of_all_generations = list(set(haplotype_to_generation.values()))
for i, generation in enumerate(list_of_all_generations):
    generation_color[generation] = generation_colormap(np.linspace(0, 1, len(list_of_all_generations)))[i]
print("Generation colors", generation_color)
for haplotype in cols:
    if haplotype in haplotype_to_generation:
        generation = haplotype_to_generation[haplotype]
        cols[haplotype] = matplotlib.colors.to_hex(generation_color[generation])


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

def add_generation_info_to_xlabels(plt, haplotype_to_generation):
    xlabels = plt.gca().get_xticklabels()
    for label in xlabels:
        haplotype = int(label.get_text())
        if haplotype in haplotype_to_generation:
            generation = int(haplotype_to_generation[haplotype])
            label.set_color(generation_color[generation])
            label.set_text(label.get_text() + " gen " + str(generation))
    plt.gca().set_xticklabels(xlabels)

add_generation_info_to_xlabels(plt, haplotype_to_generation)
ax_dend.set_ylabel('Haplotype age/generations',fontsize=24)
ax_dend.set_title('Haplotype clusters',fontsize=24)



# Plot dendrogram and colour branches
ax_dend_2 = fig.add_subplot(gs[1, 0])

dflt_col = "#808080"
different_gens_col = "#000000"

link_cols = {}
for i, i12 in enumerate(Z[:,:2].astype(int)):
    c1, c2 = (link_cols[x] if x > len(Z) else cols[x] for x in i12)

    this_col = str()
    if c1 == c2:
        this_col = c1
    elif c1 != dflt_col and c2 != dflt_col:
        this_col = different_gens_col
    else:
        this_col = dflt_col
    link_cols[i+1+len(Z)] = this_col

sns.despine(ax=ax_dend_2, offset=5, bottom=True, top=True)
dd = sch.dendrogram(Z,link_color_func=lambda x: link_cols[x], ax=ax_dend_2)

ls = []
for leaf, leaf_color in zip(plt.gca().get_xticklabels(), dd["leaves_color_list"]):
    leaf.set_color(cols[int(leaf.get_text())])
    ls.append(int(leaf.get_text()))

add_generation_info_to_xlabels(plt, haplotype_to_generation)

for i in range(len(generation_color.keys()) + 1):
    ax_dend_2.plot([]) # hack to allow us to draw legend
ax_dend_2.set_ylabel('Haplotype age/generations',fontsize=24)
ax_dend_2.legend([str(int(key)) for key in generation_color.keys()] + ['mix'], labelcolor=[generation_color[key] for key in generation_color.keys()] + ['#000000'], handlelength=0, fontsize='xx-small')

plt.savefig(output_path)