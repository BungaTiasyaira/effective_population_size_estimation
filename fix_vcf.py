# script that deletes all position == 35001 and M = 0 mutations in vcf files

import os

# getting all the directories of where we keep vcf files 
sample_list = [10,100]
AF_list = [0.5, 0.7]
VCF_directories = []

for i in range(1,11):
    for sample_size in sample_list:
        for AF in AF_list:
            directory = f"{i}/VCF/{AF}/{sample_size}"
            VCF_directories.append(directory)

output_file = "vcf_file_paths.txt"
all_vcf_paths = []

# to a text file and a list we output all the specific filepaths of the vcfs
with open(output_file, "w") as outfile:
    for directory in VCF_directories:
        if os.path.isdir(directory):
            for filename in os.listdir(directory):
                if filename.endswith(".vcf"):
                    full_path = os.path.join(directory, filename)
                    outfile.write(full_path + "\n")
                    all_vcf_paths.append(full_path)
        else:
            print(f"Warning: Directory does not exist - {directory}")

# fill in new_vcf paths
new_paths = []

for full_path in all_vcf_paths:
    parts = full_path.split(os.sep)
    # Replace 'VCF' with 'VCF_NEW'
    parts[1] = "VCF_NEW"
    new_path = os.path.join(*parts)
    new_paths.append(new_path)

# now we edit all these vcf files 

for i in range(len(all_vcf_paths)):
    # get pair of input and output vcf file paths
    input_vcf = all_vcf_paths[i]
    output_vcf = new_paths[i]

    # modify each vcf

    os.makedirs(os.path.dirname(output_vcf), exist_ok=True)

    with open(input_vcf, "r") as infile, open(output_vcf, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                # Always keep header lines
                outfile.write(line)
            else:
                fields = line.strip().split('\t')
                pos = fields[1]
                info_field = fields[7]

                # Parse INFO field into a dictionary
                info_dict = dict(item.split('=') for item in info_field.split(';') if '=' in item)

                # Check if position == 35001 and MT == 0
                if pos == "35001" and info_dict.get("MT") == "0":
                    continue  # Skip this line
                else:
                    outfile.write(line)

