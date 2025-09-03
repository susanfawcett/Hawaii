#!/bin/bash
#SBATCH --job-name=getorganelle
#SBATCH --account=fc_labordia
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=getorganelle_%j.log

# Load conda environment
module load miniconda3
conda activate getorganelle

# Base directories
INPUT_DIR="/global/scratch/users/sfawcett/GetOrganelle/ArborWGSMadieae"
OUTPUT_DIR="/global/scratch/users/sfawcett/GetOrganelle/GetOrganelleOutput"
SEED_PLASTOME="/global/scratch/users/sfawcett/WGS/Venegasia_plastome_ref/Venegasia_chloroplast_genome.fasta"

# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop over R1 files
for R1 in "$INPUT_DIR"/*_R1.fastq; do
    # Get sample base name
    SAMPLE=$(basename "$R1" _R1.fastq)
    R2="${INPUT_DIR}/${SAMPLE}_R2.fastq"
    
    # Make sample output directory
    SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"
    mkdir -p "$SAMPLE_OUT"
    
    echo "Running GetOrganelle for $SAMPLE"
    
    get_organelle_from_reads.py \
        -1 "$R1" \
        -2 "$R2" \
        -o "$SAMPLE_OUT" \
	    --continue \
        -R 15 \
        -k 21,45,65,85,105 \
        -F embplant_pt \
        -s "$SEED_PLASTOME"
done

