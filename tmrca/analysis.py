import numpy as np
import scipy.cluster.hierarchy as sch
from scipy.optimize import fsolve
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import collections
from itertools import combinations
import sys, os
from analysis_trees import get_mapping_of_haplotype_to_mutation, Node, Cluster, Haplotype

# Kit this file is just the analysis part from tmrca.py I basically copied it here
def load_data_for_analysis(Z_path, cols_path):
    Z = np.load(Z_path, allow_pickle=True)
    cols = np.load(cols_path, allow_pickle=True).item() # .item() is needed to load as dict
    return Z, cols

def get_vcf_basename(filename):
    return os.path.basename(filename)[:-4]

def get_paths(vcz_filename, threshold=None, points=None, window=None):
    basename = get_vcf_basename(vcz_filename)
    basename_components = basename.split('_')
    seed_index = int(basename_components[0])
    burnin_number = int(basename_components[1])
    population_size = int(basename_components[2])
    mutation_rate = float(basename_components[3])
    recombination_rate = float(basename_components[4])
    frequency = float(basename_components[5])
    sample_size = int(basename_components[6])

    file_general_format = '_'.join(basename_components[:-1]) # exclude sample_size from filename

    if threshold is None or points is None or window is None:
        Z_path = None
        cols_path = None
        output_path = None
    else:
        tmrca_folder = f'{seed_index}/TMRCA/{frequency}/{sample_size}'
        tmrca_output_format = f'{tmrca_folder}/{file_general_format}_{sample_size}_{threshold}_{points}_{window}'
        Z_path = f'{tmrca_output_format}_Z.npy'
        cols_path = f'{tmrca_output_format}_cols.npy'

        output_folder = f'{seed_index}/TMRCA_dendograms/{frequency}/{sample_size}'
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        output_path = f'{output_folder}/{file_general_format}_{sample_size}_{threshold}_{points}_{window}_analysis.pdf'

    trees_folder = f'{seed_index}/TREES'
    trees_path = f'{trees_folder}/{file_general_format}.trees'

    sample_ID_folder = f'{seed_index}/sample_ID/{frequency}/{sample_size}'
    sample_ID_path = f'{sample_ID_folder}/{file_general_format}_{sample_size}.txt'

    return Z_path, cols_path, output_path, trees_path, sample_ID_path, population_size, mutation_rate

def run_tmrca_analysis(Z, cols, trees_path, sample_ID_path, output_path, population_size, mutation_rate, meta_path, name):
    print('Z', Z.shape, Z)
    print('cols', cols)

    haplotype_to_mutation = get_mapping_of_haplotype_to_mutation(trees_path,
                                                                    sample_ID_path)
    mutations = sorted(list(set(haplotype_to_mutation.values())), key=lambda mut: mut.get_generation())
    print('Mutations:', mutations)
    for i, mutation in enumerate(mutations):
        mutation.set_color(matplotlib.colors.to_hex(matplotlib.colormaps['gist_rainbow'](np.linspace(0, 1, len(mutations)))[i]))

    haplotypes = [Haplotype(haplotype, haplotype_to_mutation.get(haplotype)) for haplotype in cols]
    print("Haplotypes:", list(filter(lambda x: x.is_mutant(), haplotypes)))


    def add_mutation_info_to_xlabels(plt, haplotype_to_mutation):
        xlabels = plt.gca().get_xticklabels()
        for label in xlabels:
            haplotype = int(label.get_text())
            if haplotype in haplotype_to_mutation:
                mutation = haplotype_to_mutation[haplotype]
                label.set_color(mutation.get_color())
                label.set_text(label.get_text() + " gen " + str(mutation.get_generation()))
        plt.gca().set_xticklabels(xlabels)


    matplotlib.rcParams.update({'font.size': 24})
    fig = plt.figure(figsize=(30, 12))
    gs = matplotlib.gridspec.GridSpec(2, 1, hspace=0.1, wspace=1, height_ratios=(1.5,1)) # type: ignore



    # Plot dendrogram and colour branches
    ax_dend = fig.add_subplot(gs[0, 0])
    ax_dend.set_title('Haplotype clusters',fontsize=24)

    dflt_col = "#808080"
    different_gens_col = "#000000"

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
        elif c1 != dflt_col and c2 != dflt_col:
            this_col = different_gens_col
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

    oldest_node = max([node.get_age() for node in node_map.values()])
    print(f"Oldest node is {oldest_node}")

    def score(mutation):
        clusters_with_this_mutation = [cluster for cluster in clusters if mutation in cluster.get_haplotypes_by_mutation()]


        def score_no_of_clusters():
            return 1 / len(clusters_with_this_mutation)

        def score_infiltration():
            correct = 0
            total = 0
            for cluster in clusters_with_this_mutation:
                cluster_haplotypes = cluster.get_haplotypes()

                correct += len([haplotype for haplotype in cluster_haplotypes if haplotype.mutation == mutation])
                total += len(cluster_haplotypes)

            return correct / total

        def score_mrca():
            combs = [(clusters_with_this_mutation[l_idx], clusters_with_this_mutation[r_idx]) for (l_idx, r_idx) in list(combinations(range(len(clusters_with_this_mutation)), 2))]
            if not combs:
                return 1.0
            score = 0
            for cluster1, cluster2 in combs:
                mrca_node = cluster1.get_mrca_to_node(cluster2.root)
                this_score = (1 - mrca_node.get_age() / oldest_node) ** 2
                score += this_score
                #print(cluster1, cluster2, mrca_node.get_age(), this_score)

            return score / len(combs)

        no_of_clusters_score = score_no_of_clusters()
        infiltration_score = score_infiltration()
        mrca_score = score_mrca()

        print("-" * 60)
        print(f"Scoring {mutation}")
        print(f"Relevant clusters: {clusters_with_this_mutation}")
        print("-" * 30)
        print(f"No of clusters score: {no_of_clusters_score:.3}")
        print(f"Infiltration score {infiltration_score:.3}")
        print(f"MRCA score {mrca_score:.3}")
        print("-" * 30)

        total_score = 0.70 * no_of_clusters_score + 0.15 * infiltration_score + 0.15 * mrca_score
        print(f"Total score {total_score:.3}")

        return no_of_clusters_score, infiltration_score, mrca_score, total_score


    print(clusters)
    for cluster in clusters:
        print(cluster.get_haplotypes_by_mutation())

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

    ax_dend_2 = fig.add_subplot(gs[1, 0])
    score_matrix = []
    analysis_no_of_clusters_score, analysis_infiltration_score, analysis_mrca_score, analysis_total_score = 0, 0, 0, 0
    analysis_num_of_haplotypes = sum([len(mutation.haplotypes) for mutation in mutations])
    for mutation in mutations:  
        no_of_clusters_score, infiltration_score, mrca_score, total_score = score(mutation)
        no_of_haplotypes = len(mutation.haplotypes)
        weight = no_of_haplotypes / analysis_num_of_haplotypes

        analysis_no_of_clusters_score += no_of_clusters_score * weight
        analysis_infiltration_score += infiltration_score * weight
        analysis_mrca_score += mrca_score * weight
        analysis_total_score += total_score * weight

        score_matrix.append([int(mutation.get_generation()), no_of_haplotypes, round(no_of_clusters_score, 2), round(infiltration_score, 2), round(mrca_score, 2), f'{round(no_of_haplotypes * total_score, 2)}/{float(no_of_haplotypes)}'])
    
    ax_dend_2.table(score_matrix, colLabels=['Mutation generation', 'Haplotype count', 'No of clusters score', 'Infiltration score', 'Mrca score', 'Final score'], bbox=[0, 0.1, 1, 0.9], cellColours=[[f'{mutation.get_color()}88']*len(score_matrix[-1]) for mutation in mutations], fontsize=13)
    ax_dend_2.text(0.5, 0, f'{analysis_total_score}', ha='center')
    ax_dend_2.axis('off')


    # print('TMRCA dendogram output to', output_path)
    # plt.savefig(output_path)

    # Estimate num of independent origins
    def equation(theta, num_independent_origins, num_mutant_haplotypes):
        return theta * np.log(1 + num_mutant_haplotypes / theta) - num_independent_origins

    def solve_theta(num_independent_origins, num_mutant_haplotypes):
        initial_guess = 1.0
        theta_solution, = fsolve(equation, initial_guess, args=(num_independent_origins, num_mutant_haplotypes))
        return theta_solution

    def calculate_Ne_from_n(num_independent_origins, num_mutant_haplotypes, mutation_rate):
        theta = solve_theta(num_independent_origins, num_mutant_haplotypes)
        Ne = theta/(4*mutation_rate)
        return Ne

    population_haplotype_to_mutation = get_mapping_of_haplotype_to_mutation(trees_path, None)
    population_haplotypes = [Haplotype(haplotype, population_haplotype_to_mutation.get(haplotype)) for haplotype in range(0, population_size*2)]
    n_d_0 = len(list(filter(lambda c: len(c.get_haplotypes()) >= 0,clusters)))
    n_d_1 = len(list(filter(lambda c: len(c.get_haplotypes()) >= 1,clusters)))
    n_d_2 = len(list(filter(lambda c: len(c.get_haplotypes()) >= 2,clusters)))
    n_d_3 = len(list(filter(lambda c: len(c.get_haplotypes()) >= 3,clusters)))
    n_d_4 = len(list(filter(lambda c: len(c.get_haplotypes()) >= 4,clusters)))

    n_s = len(mutations)
    n_p = len(set(population_haplotype_to_mutation.values()))

    Ne_d_0 = calculate_Ne_from_n(n_d_0, len(list(filter(lambda h: h.is_mutant(), haplotypes))), mutation_rate) 
    Ne_d_1 = calculate_Ne_from_n(n_d_1, len(list(filter(lambda h: h.is_mutant(), haplotypes))), mutation_rate)
    Ne_d_2 = calculate_Ne_from_n(n_d_2, len(list(filter(lambda h: h.is_mutant(), haplotypes))), mutation_rate)
    Ne_d_3 = calculate_Ne_from_n(n_d_3, len(list(filter(lambda h: h.is_mutant(), haplotypes))), mutation_rate)
    Ne_d_4 = calculate_Ne_from_n(n_d_4, len(list(filter(lambda h: h.is_mutant(), haplotypes))), mutation_rate)

    Ne_s = calculate_Ne_from_n(n_s, len(list(filter(lambda h: h.is_mutant(), haplotypes))), mutation_rate)
    Ne_p = calculate_Ne_from_n(n_p, len(list(filter(lambda h: h.is_mutant(), population_haplotypes))), mutation_rate)
    Ne_actual = population_size

    print(n_d_0, n_d_1, n_d_2, n_d_3, n_d_4, n_s, n_p, Ne_d_0, Ne_d_1, Ne_d_2, Ne_d_3, Ne_d_4, Ne_s, Ne_p, Ne_actual)

    with open(meta_path, "a") as file:
        # Lock the meta file
        os.lockf(file.fileno(), os.F_LOCK, 0)
        # file.write(f'{name},{threshold},{points},{window},{analysis_no_of_clusters_score},{analysis_infiltration_score},{analysis_mrca_score},{analysis_total_score},{n_d},{n_s},{n_p},{Ne_d},{Ne_s},{Ne_p},{Ne_actual},{num_mutant_haplotypes}\n')
        file.write(f'{name},{analysis_no_of_clusters_score},{analysis_infiltration_score},{analysis_mrca_score},{analysis_total_score},{n_d_0},{n_d_1},{n_d_2},{n_d_3},{n_d_4},{n_s},{n_p},{Ne_d_0},{Ne_d_1},{Ne_d_2},{Ne_d_3},{Ne_d_4},{Ne_s},{Ne_p},{Ne_actual},{len(list(filter(lambda h: h.is_mutant(), haplotypes)))},{len(list(filter(lambda h: h.is_mutant(), population_haplotypes)))} \n')
        # Release the lock
        os.lockf(file.fileno(), os.F_ULOCK, 0)

if __name__ == "__main__":
    try:
        vcz_filename = sys.argv[1]
        threshold = float(sys.argv[2])
        points = int(sys.argv[3])
        window = int(sys.argv[4])
    except IndexError:
        print("Usage: python3 analysis.py <vcz_filename> <threshold> <points> <window>")
        exit()

    Z_path, cols_path, output_path, trees_path, sample_ID_path, population_size, mutation_rate = get_paths(vcz_filename, threshold, points, window)

    Z, cols = load_data_for_analysis(Z_path, cols_path)

    run_tmrca_analysis(Z, cols, trees_path, sample_ID_path, output_path, population_size, mutation_rate, "meta_Ne.csv", get_vcf_basename(vcz_filename))