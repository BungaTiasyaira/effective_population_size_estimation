#!/bin/bash

#PBS -N TMRCA
#PBS -j oe
#PBS -J 1-19

#PBS -m ae

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=4:mem=15gb

## NB values for ncpus and mem are allocated
## to each node (specified by select=N)
##
## to direct output to cwd, use $PBS_O_WORKDIR:
## specify LOGFILE found in ~/ during execution then moved to cwd on job completion
##
cd $PBS_O_WORKDIR
JOBNUM=`echo $PBS_JOBID | sed 's/\..*//'`
LOGFILE=${PBS_JOBNAME}.o${JOBNUM}.i${PBS_ARRAY_INDEX}

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

echo "PBS job array index is $PBS_ARRAY_INDEX"

FILE="active_zars.txt"
LINE_NUM="$PBS_ARRAY_INDEX"

LINE_NUM=$((LINE_NUM+1)) # skip column names

if [[ ! -f "$FILE" ]]; then
    echo "Error: File '$FILE' does not exist."
    exit 1
fi

if ! [[ "$LINE_NUM" =~ ^[0-9]+$ ]]; then
    echo "Error: Line number must be a positive integer."
    exit 1
fi

LINE=$(sed -n "${LINE_NUM}p" "$FILE")

if [[ -z "$LINE" ]]; then
    echo "Error: Line $LINE_NUM does not exist in '$FILE'."
    exit 1
fi


IFS=',' read -ra words <<< "$LINE"
echo ${words[0]} ${words[1]} ${words[2]} ${words[3]}

python3 tmrca/tmrca.py ${words[0]} ${words[1]} ${words[2]} ${words[3]}

echo "Finished running TMRCA"