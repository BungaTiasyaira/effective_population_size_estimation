# import modules 
import tskit
import numpy as np
import pandas as pd
from collections import Counter
import re


# function to read parameters file and get file names 
def parse_file(file_path):
    N = []
    mut_rate = []
    rec_rate = []
    burnin_num = []
    
    try:
        with open(file_path, 'r') as file:
            # Skip the header
            next(file)
            
            for line in file:
                parts = line.strip().split(',')
                if len(parts) == 4:
                    N.append(int(parts[0]))
                    mut_rate.append(float(parts[1]))
                    rec_rate.append(float(parts[2]))
                    burnin_num.append(int(parts[3]))
    
        return N, mut_rate, rec_rate, burnin_num
    
    except FileNotFoundError:
        print("File not found. Please check the file path.")
        return [], [], [], []
    except ValueError as e:
        print(f"Data format error: {e}")
        return [], [], [], []


# function to read sample ID file and convert to an array of integers
def file_to_array(file_path):
    try:
        with open(file_path, 'r') as file:
            # Read lines, strip spaces, and convert to integers
            array = [int(line.strip()) for line in file if line.strip().isdigit()]
        return array
    except FileNotFoundError:
        print("File not found. Please check the file path.")
        return []
    except ValueError:
        print("File contains non-numeric data.")
        return []


# function that fills dictionrary
def get_pair(allele_frequency, sample_size, dictionary):
    SampleID_filename = f"{Seed_ID}/sample_ID/{allele_frequency}/{sample_size}/{Seed_ID}_{burnin_num[i]}_{N[i]}_{mut_rate[i]}_{rec_rate[i]}_{allele_frequency}_{sample_size}.txt"
    tree_filename = f"{Seed_ID}/TREES/{Seed_ID}_{burnin_num[i]}_{N[i]}_{mut_rate[i]}_{rec_rate[i]}_{allele_frequency}.trees"
    dictionary[SampleID_filename] = tree_filename

def extract_values(sampleID_filepath):
    pattern = r"(\d+)/sample_ID/\d+\.\d+/(\d+)/(\d+)_(\d+)_(\d+)_([0-9]+\.[0-9]+)"

    # Iterate over the file paths and extract values using regex
    match = re.search(pattern, sampleID_filepath)
    if match:
        # Extract values from the regex match groups
        seed_id = match.group(1)
        sample_size = match.group(2)  # Now captures the correct sample size
        burnin_number = match.group(4)
        mutation_rate = match.group(6)  # Mutation rate is the 6th group
        return seed_id, burnin_number, sample_size, mutation_rate
    else:
        return None

# function that does the analysis and puts file in output in a set directory 
def analyse_sampleID(sampleID_filepath, dictionary):

    # use function above to access sampleID file and turn into array
    sample_IDs = file_to_array(sampleID_filepath)

    # Try to load the .trees file using tskit
    try:
        ts = tskit.load(dictionary[sampleID_filepath])
    except FileNotFoundError:
        print(f"File not found: {dictionary.get(sampleID_filepath, 'Unknown')}. Skipping to the next file.")
        return
    except Exception as e:
        print(f"An error occurred while loading the .trees file for {sampleID_filepath}: {e}")
        return

    # initialize empty DataFrame with column names-- to store data as csv
    # we will add to this later on by using df.loc[x] = ['Genome ID', 'Genotype ID', 'Allele']
    df = pd.DataFrame(columns=['Genome ID', 'Genotype ID', 'Allele'])

    sample_genomes = []
    sample_genotypes = []

    for individual in ts.individuals():
        if individual.id in sample_IDs:
            sample_genomes.append(int(individual.nodes[0]))
            sample_genomes.append(int(individual.nodes[1]))

    # Based on our sampled genomes, get the genotypes for said genome and what allele corresponds to it!
    for variant in ts.variants():
        if variant.site.position == 35000:
            # Get the indices of sample_genomes elements in variant_samples (should be the same but just to be sure)
            indices = [np.where(variant.samples == num)[0][0] for num in sample_genomes]
            for indice in indices:
                df.loc[len(df)] = [indice, variant.genotypes[indice], variant.alleles[variant.genotypes[indice]]] # filling in our dataframe
                sample_genotypes.append(variant.genotypes[indice]) # filling in our genotypes list 

    # Now I want to count the number of lineages present in the sample, along with their frequencies. Best to output this as another .csv file
    df_frequencies = pd.DataFrame(columns=['Genotype ID', 'Frequency (%)']) 

    unique_sample_genotypes = list(set(sample_genotypes)) # only genotypes present in sample
    count = Counter(sample_genotypes)

    for variant in ts.variants():
        if variant.site.position == 35000:
            for i in range (len(unique_sample_genotypes)):
                number = count[unique_sample_genotypes[i]]
                frequency = (count[unique_sample_genotypes[i]]/len(sample_genotypes))*100
                df_frequencies.loc[len(df_frequencies)] = [unique_sample_genotypes[i], frequency]
    
    # adding to the all csv dataframes

    if '0.5' in sampleID_filepath:
        result = extract_values(sampleID_filepath)
        seed_id = result[0]
        burnin_number = result[1]
        sample_size = result[2]
        Independent_origins = len([x for x in unique_sample_genotypes if x != 0])
        Mutants_in_Sample = len([y for y in sample_genotypes if y != 0])
        Mutation_rate = result[3]
        df_all_05.loc[len(df_all_05)] = [seed_id, burnin_number, sample_size, Independent_origins, Mutants_in_Sample, Mutation_rate]  # Add a row directly

    if '0.7' in sampleID_filepath:
        result = extract_values(sampleID_filepath)
        seed_id = result[0]
        burnin_number = result[1]
        sample_size = result[2]
        Independent_origins = len([x for x in unique_sample_genotypes if x != 0])
        Mutants_in_Sample = len([y for y in sample_genotypes if y != 0])
        Mutation_rate = result[3]
        df_all_07.loc[len(df_all_07)] = [seed_id, burnin_number, sample_size, Independent_origins, Mutants_in_Sample, Mutation_rate]  # Add a row directly

    # output to csv
    last_dot_index = sampleID_filepath.rfind('.')

    df.to_csv(f"{sampleID_filepath[0:last_dot_index]}_Alleles.csv", index=False) 
    df_frequencies.to_csv(f"{sampleID_filepath[0:last_dot_index]}_Frequencies.csv", index=False) 
    print(f"CSV files for {sampleID_filepath} created successfully!")

# access parameters file and return 4 arrays 
param_file_path = 'parameter_combinations.txt'  
N, mut_rate, rec_rate, burnin_num = parse_file(param_file_path)

print("N:", N)
print("Mutation Rates:", mut_rate)
print("Recombination Rates:", rec_rate)
print("Burnins:", burnin_num)

# create seed array from 1-10
Seed_IDs = [1,2,3,4,5,6,7,8,9,10]


# initialize empty DataFrames with column names-- to store data as csv
df_all_05 = pd.DataFrame(columns=['Seed', 'Burnin No.', 'Sample Size', 'Number of Independent Origins', 'Mutants in Sample', 'Mutation_rate'])
df_all_07 = pd.DataFrame(columns=['Seed', 'Burnin No.', 'Sample Size', 'Number of Independent Origins', 'Mutants in Sample', 'Mutation_rate'])


# Dictionary containing file pairs (sample_IDs & .tree file)
# Key = sample_ID filepath
# Value = .trees filepath 
file_pairs_05 = {}
file_pairs_07 = {}

# write a loop that goes over the arrays and gets all sample ID & trees filename pairs (in form of dictionary!)
for Seed_ID in Seed_IDs:
    for i in range(len(N)):
        if N[i] >= 10000:

            get_pair(0.5, 1000, file_pairs_05)
            get_pair(0.7, 1000, file_pairs_07)

        if N[i] >= 1000:
            
            get_pair(0.5, 100, file_pairs_05)
            get_pair(0.7, 100, file_pairs_07)

        if N[i] >= 100:

            get_pair(0.5, 10, file_pairs_05)
            get_pair(0.7, 10, file_pairs_07)

# run analyse_sampleID function for both dictionaries (0.5 and 0.7)
# the csv files will be in the same directory as where the Sample IDs are with the same file name (except with Alleles/Frequencies at the end)

for key in file_pairs_05:
    analyse_sampleID(key, file_pairs_05)

for key in file_pairs_07:
    analyse_sampleID(key, file_pairs_07)

# save csv all

df_all_05.to_csv("IndependentOrigins_05.csv", index=False)
df_all_07.to_csv("IndependentOrigins_07.csv", index=False)