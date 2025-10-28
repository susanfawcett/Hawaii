#!/bin/bash

#SBATCH --job-name=hybassemble
#SBATCH --account=fc_labordia
#SBATCH --partition=savio3
#SBATCH --cpus-per-task=4
#SBATCH --time=68:00:00
#SBATCH --output=/global/scratch/users/sfawcett/HybPiper/logs/hybassemble_%j.out
#SBATCH --error=/global/scratch/users/sfawcett/HybPiper/logs/hybassemble_%j.err

date

module load anaconda3
source /global/home/users/sfawcett/miniconda3/etc/profile.d/conda.sh
conda activate hybpiper

WDIR=/global/scratch/users/sfawcett/HybPiper/Loulu_Data/trimmed_reads
TARGETS=/global/scratch/users/sfawcett/HybPiper/mega353_fixed.fasta
NAMELIST=/global/scratch/users/sfawcett/HybPiper/namelist.txt

cd $WDIR
mkdir -p /global/scratch/users/sfawcett/HybPiper/logs

while read name; do
    echo "Running HybPiper on $name..."
    hybpiper assemble -t_dna $TARGETS -r ${name}*.fastq --prefix $name --bwa
    echo "Finished $name"
done < $NAMELIST

conda deactivate

echo "ALL DONE"
date
