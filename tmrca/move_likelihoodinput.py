# moving output from MEOWSS meta.csv to get_likelihood_input.csv

import csv

# Input and output file names
input_file = 'metas/MEOWSS_meta.csv'
output_file = 'tmrca/get_likelihood_input.csv'

# Open the input file
with open(input_file, newline='') as infile:
    reader = csv.DictReader(infile)

    # Prepare output data
    output_rows = []

    for row in reader:
        name = row['name']
        n1 = row['n_d_1'] 
        n2 = row["n_d_2"] 
        n3 = row["n_d_3"] 
        n4 = row["n_d_4"]           
        Ns = row['mut_haplotypes_sample']  
        label = name

        if '16' in row['name'] or '17' in row['name']:
            mu = '0.0025'    
        if '13' in row['name'] or '14' in row['name'] or '15' in row['name']:
            mu = '0.00025' 

        output_rows.append({
            'n1': n1,
            'n2': n2,
            'n3': n3,
            'n4': n4,             
            'Ns': Ns,
            'mu': mu,
            'label': label,
        })

# Write to the output file
with open(output_file, 'w', newline='') as outfile:
    fieldnames = ['n1','n2','n3','n4', 'Ns', 'mu', 'label']
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)

    writer.writeheader()
    for row in output_rows:
        writer.writerow(row)

print(f"Data written to {output_file}")
