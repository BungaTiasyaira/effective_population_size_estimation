#!/bin/bash

#PBS -N phasedibd_gettingruntimes
#PBS -j oe
#PBS -J 1-15

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
START_TIME=$(date +%s)

# Read  line from the list of vcf paths
VCF_PATH=$(sed -n "${PBS_ARRAY_INDEX}p" get_vcfs_for_SINGER.txt)

# Run the Python script with the selected vcf path
echo "Now processing: $VCF_PATH"
python3 tmrca/tmrca_phaseibd.py "$VCF_PATH.vcf"
echo "finished making Z and and dendogram for $VCF_PATH using phaseibd"

TMRCA_END_TIME=$(date +%s)
TMRCA_ELAPSED_TIME=$((TMRCA_END_TIME - START_TIME))
echo "Total runtime tmrca: $(date -ud "@$TMRCA_ELAPSED_TIME" +'%H hours %M minutes %S seconds')"

