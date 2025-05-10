import os

def parse_log_path(log_path):
    basename = os.path.basename(log_path)[:-8]
    print(basename)
    basename_components = basename.split('_')
    seed_index = int(basename_components[0])
    burnin_number = int(basename_components[1])
    population_size = int(basename_components[2])
    mutation_rate = float(basename_components[3])
    recombination_rate = float(basename_components[4])

    file_general_format = '_'.join(basename_components[:-1]) # exclude sample_size from filename

    return seed_index, burnin_number, population_size, mutation_rate, recombination_rate

def parse_vcz_path(vcz_path):
    basename = os.path.basename(vcz_path)[:-4]
    basename_components = basename.split('_')
    seed_index = int(basename_components[0])
    burnin_number = int(basename_components[1])
    population_size = int(basename_components[2])
    mutation_rate = float(basename_components[3])
    recombination_rate = float(basename_components[4])
    frequency = float(basename_components[5])
    sample_size = int(basename_components[6])

    file_general_format = '_'.join(basename_components[:-1]) # exclude sample_size from filename

    return seed_index, burnin_number, population_size, mutation_rate, recombination_rate, frequency, sample_size, file_general_format

def get_tmrca_paths(basename, threshold, points, window, makeDirs=True):
    seed_index, burnin_number, population_size, mutation_rate, recombination_rate, frequency, sample_size, file_general_format = parse_full_path(basename)

    tmrca_folder = f'{seed_index}/TMRCA/{frequency}/{sample_size}'
    tmrca_output_format = f'{tmrca_folder}/{file_general_format}_{sample_size}_{threshold}_{points}_{window}'
    Z_path = f'{tmrca_output_format}_Z.npy'
    cols_path = f'{tmrca_output_format}_cols.npy'

    trees_folder = f'{seed_index}/TREES'
    trees_path = f'{trees_folder}/{file_general_format}.trees'

    sample_ID_folder = f'{seed_index}/sample_ID/{frequency}/{sample_size}'
    sample_ID_path = f'{sample_ID_folder}/{file_general_format}_{sample_size}.txt'

    output_folder = f'{seed_index}/TMRCA_dendograms/{frequency}/{sample_size}'
    if makeDirs and not os.path.exists(output_folder):
        os.makedirs(output_folder)
    output_path = f'{output_folder}/{file_general_format}_{sample_size}_{threshold}_{points}_{window}_analysis.pdf'

    return Z_path, cols_path, trees_path, sample_ID_path, output_path