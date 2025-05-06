#!/bin/bash

# file ba1hap.sh
#

#SBATCH --job-name=hap_1
#SBATCH --account=fc_labordia
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=71:58:59

## commands:
echo
echo "START"; date
echo
echo "conda: activate sorter"

 # Initialize conda
source /global/home/users/sfawcett/miniconda3/etc/profile.d/conda.sh
conda activate /global/home/users/sfawcett/miniconda3/envs/SORTER2

export PYTHON=/global/home/users/sfawcett/miniconda3/envs/SORTER2/bin/python



export WDIR=/global/scratch/users/sfawcett/SORTER2/Pritchardia/Raw/SORTER2_Pritchardia/
cd $WDIR
echo
echo "WDIR = $PWD"
echo

echo $SLURM_JOB_NODELIST | sed s/\,/\\n/g > hostfile

CMD="/global/home/users/sfawcett/miniconda3/envs/SORTER2/bin/python /global/scratch/users/akostic/sorter/bin/hap1.py -wd $WDIR -cpref /global/scratch/users/sfawcett/Data/Pritchardia/Plastomes/Pritchardia_pacifica_plastome.fasta -c 80 -d 10"

# print the command line
echo "batch:CMD=$CMD"

# run it now
$CMD

echo
date
echo "DONE"
