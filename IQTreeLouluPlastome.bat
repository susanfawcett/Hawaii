#!/bin/bash

# file MafftPritchardiaPlastome.sh

#SBATCH --job-name=IQTreePritchardiaPlastome
#SBATCH --account=fc_labordia
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=71:58:59

# Directory containing the aligned files - first untrimmed, second trimal -strictplus after taxa with <50% coverage removed
#ALIGNMENT_FILE="/global/scratch/users/sfawcett/Haplominer/LouluPlastomeAlignment/LouluPlastomes_4May2025.fasta"
ALIGNMENT_FILE="/global/scratch/users/sfawcett/Haplominer/LouluPlastomeAlignment/trimmed_strictplus.fasta"

# Output directory for IQ-TREE result
IQTREE_OUTPUT_DIR="/global/scratch/users/sfawcett/Haplominer/IQTreeLouluPlastome/LouluPlastome_strictplus"

# Create the output directory if it doesn't exist
mkdir -p $IQTREE_OUTPUT_DIR

# Number of bootstrap replicates
BOOTSTRAPS=1000

# Number of threads to use
THREADS=AUTO

# Move to the output directory
cd $IQTREE_OUTPUT_DIR
  
# Run IQ-TREE on the alignment file
  iqtree3 -s $ALIGNMENT_FILE -m TEST -bb $BOOTSTRAPS -T $THREADS -safe --prefix louluplastome

echo
echo "IQ-TREE analyses completed!"
date
