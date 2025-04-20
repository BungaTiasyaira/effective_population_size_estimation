#!/bin/bash

#PBS -N syaii_simulations
#PBS -j oe
#PBS -k oe

#PBS -m ae

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=24:mem=15gb

## NB values for ncpus and mem are allocated
## to each node (specified by select=N)
##
## to direct output to cwd, use $PBS_O_WORKDIR:
## specify LOGFILE found in ~/ during execution then moved to cwd on job completion
##
cd $PBS_O_WORKDIR
JOBNUM=`echo $PBS_JOBID | sed 's/\..*//'`
LOGFILE=${PBS_JOBNAME}.o${JOBNUM}

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

filename="parameter_combinations.txt"
seed_filename="seeds.txt"

# Check if file exists
if [ ! -f "$filename" ]; then
    echo "Parameter File not found!"
    exit 1
fi

if [ ! -f "$seed_filename" ]; then
    echo "Seed File not found!"
    exit 1
fi

# Read file line by line
while IFS= read -r line; do
    IFS=',' read -ra words <<< "$line"
    seed_id=${words[0]}
    seed=${words[1]}
    tail -n +2 "$filename" | while IFS= read -r line; do
        # Split by comma and read each word
        IFS=',' read -ra words <<< "$line"
        echo Running simulation with population_size=${words[0]} mut_rate=${words[1]} recomb_rate=${words[2]} burnin_number=${words[3]} seed=${seed}
        slim -d population_size=${words[0]} -d mut_rate=${words[1]} -d recomb_rate=${words[2]} -d burnin_number=${words[3]} -d seed=${seed} -d seed_id=${seed_id} SS_all.slim
        echo "---" # Separator for clarity
    done
done < "$seed_filename"

## move LOGFILE to cwd
##
echo "Finished running simulations"
