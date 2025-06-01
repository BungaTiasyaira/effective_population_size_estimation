#!/bin/bash

#PBS -N vcf2zarr
#PBS -j oe
#PBS -J 1-10

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

rm -rf "$PBS_ARRAY_INDEX/ZARR"
cp -r "$PBS_ARRAY_INDEX/VCF" "$PBS_ARRAY_INDEX/ZARR"
DIR="$PBS_ARRAY_INDEX/ZARR"

find "$DIR" -type f | while read -r file; do
    echo "$file"
    bgzip "$file"
    tabix "$file.gz"
    python3 -m bio2zarr vcf2zarr explode "$file.gz" "${file%????}.icf"
    python3 -m bio2zarr vcf2zarr encode "${file%????}.icf" "${file%????}.vcz"
    rm "$file.gz"
    rm "$file.gz.tbi"
    rm -rf "${file%????}.icf"
done

echo "Finished running vcf2zarr"
