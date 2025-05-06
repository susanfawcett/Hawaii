#!/bin/bash

# file batch_sorter_1A.sh
#

#SBATCH --job-name=sorter1A
#SBATCH --account=fc_labordia
#SBATCH --partition=savio3_bigmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=71:58:59

export SCR="/global/scratch/users/sfawcett"


# Initialize conda
source /global/home/users/sfawcett/miniconda3/etc/profile.d/conda.sh
conda activate SORTER2

## commands:
echo
echo "START"; date
echo
echo "conda: activate SORTER2"

export WDIR=/global/scratch/users/sfawcett/SORTER2/Pritchardia/Raw
cd $WDIR
echo
echo "WDIR = $PWD"
echo

echo $SLURM_JOB_NODELIST |sed s/\,/\\n/g > hostfile

export JOBS_PER_NODE=$(( $SLURM_CPUS_ON_NODE / $SLURM_CPUS_PER_TASK ))

echo "SLURM_CPUS_ON_NODE = $SLURM_CPUS_ON_NODE"
echo "SLURM_CPUS_PER_TASK = $SLURM_CPUS_PER_TASK"
echo "JOBS_PER_NODE = $JOBS_PER_NODE"
echo

python /global/scratch/users/akostic/sorter/bin/1a.py -n Pritchardia -trim T -spades T
 
# print the command line
echo "batch:CMD=$CMD"

# run it now
$CMD

echo
date
echo "DONE"
