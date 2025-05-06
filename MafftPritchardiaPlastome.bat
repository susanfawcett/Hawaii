#!/bin/bash

# file MafftPritchardiaPlastome.sh

#SBATCH --job-name=mafftPritchardiaPlastome
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

# Directory containing the plastome file
ASSEMBLY="/global/scratch/users/sfawcett/Haplominer/PritchardiaPlastomesRenamed/all_loulu_plastomes_unaligned.fasta"

# Reference plastome file
REFERENCE="/global/scratch/users/sfawcett/Data/Pritchardia/Plastomes/Pritchardia_pacifica_plastome.fasta"

# Output directory for aligned files
OUTPUT_DIR="/global/scratch/users/sfawcett/Haplominer/LouluPlastomeAlignment"

# Create the output directory if it doesn't exist
mkdir -p $OUTPUT_DIR
  
# Perform MAFFT reference-guided alignment using --add
mafft --thread 40 --auto --keeplength --add "$ASSEMBLY" "$REFERENCE" > "$OUTPUT_DIR/aligned_plastomes.fasta"


echo "Alignments completed!"
