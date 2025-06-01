#!/bin/bash

#PBS -N singer_mantku
#PBS -j oe
#PBS -J 1-100

#PBS -m ae

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=4:mem=15gb

## NB values for ncpus and mem are allocated
## to each node (specified by select=N)
##
## to direct output to cwd, use $PBS_O_WORKDIR:
## specify LOGFILE found in ~/ during execution then moved to cwd on job completion
##
cd $PBS_O_WORKDIR
JOBNUM=`echo $PBS_JOBID | sed 's/\..*//'`
LOGFILE=${PBS_JOBNAME}.o${JOBNUM}_${PBS_ARRAY_INDEX}

#########################################
##                                     ##
## Output some useful job information. ##
##                                     ##
#########################################

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: job number is $JOBNUM
echo PBS: logfile is $LOGFILE
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
#echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

## load common modules as standard
##
module load anaconda3/personal vcftools
source activate slim_fyp

# get vcf file path
VCF_PATH=$(sed -n "$((PBS_ARRAY_INDEX))p" get_vcfs_for_SINGER.txt)
VCF_FILENAME=$(basename "$VCF_PATH")

# Extract parameters from the VCF filename
# Format: 5/VCF_NEW/0.7/100/5_13_1000_0.00025_2.5e-05_0.7_100
IFS='_' read -r seed burnin N m r af samplesize_vcf <<< "${VCF_FILENAME%.vcf}"
ratio=$(python3 -c "print(float($r) / float($m))")

OUTPUT_DIR="singer_outputs/${VCF_FILENAME%.vcf}"
mkdir -p "$OUTPUT_DIR/outputs"
mkdir -p "$OUTPUT_DIR/trees"

./singer/singer_master \
    -m $m \
    -ratio $ratio \
    -vcf $VCF_PATH \
    -output "$OUTPUT_DIR/outputs/output" \
    -start 0 \
    -end 70001

echo "Finished running singer for seed: $seed, burnin: $burnin, AF: $af, mut rate: $m, rec rate: $r"
echo "converting to trees"

./singer/convert_to_tskit -input "$OUTPUT_DIR/outputs/output" -output "$OUTPUT_DIR/trees/tree" -start 0 -end 70001 -step 1

# echo "finished converting to trees"
echo "Making dendograms"

for i in $(seq 0 99); do
    FILE="$OUTPUT_DIR/trees/tree_${i}.trees"
    python tmrca/understandingSINGER_output.py "$FILE" "$i"
    echo "finished scoring $FILE"
done 

