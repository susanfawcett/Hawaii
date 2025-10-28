#!/bin/bash

#SBATCH --job-name=checktargets
#SBATCH --account=fc_labordia
#SBATCH --partition=savio3
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=/global/scratch/users/sfawcett/HybPiper/logs/targetcheck_%j.out
#SBATCH --error=/global/scratch/users/sfawcett/HybPiper/logs/targetcheck_%j.err

module load anaconda3
source /global/home/users/sfawcett/miniconda3/etc/profile.d/conda.sh
conda activate hybpiper

hybpiper check_targetfile -t_dna mega353.fasta
