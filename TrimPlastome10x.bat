#!/bin/bash
#SBATCH --job-name=TrimPlastomes10x
#SBATCH --account=fc_labordia
#SBATCH --partition=savio3
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/global/scratch/users/sfawcett/SORTER2/Pritchardia/Plastome10x/logs/Trim10x_%j.out
#SBATCH --error=/global/scratch/users/sfawcett/SORTER2/Pritchardia/Plastome10x/logs/Trim10x_%j.err

echo "START TRIMMING"
date

# --- Environment setup ---
source /global/home/users/sfawcett/miniconda3/etc/profile.d/conda.sh
conda activate SORTER2

# --- Directories ---
BASE_DIR=/global/scratch/users/sfawcett/SORTER2/Pritchardia
ASSEMBLY_DIR=$BASE_DIR/Raw/SORTER2_Pritchardia
PLASTOME_DIR=$BASE_DIR/Plastome10x/LouluPlastomes
OUT_DIR=$BASE_DIR/Plastome10x/Trimmed10x
LOG_DIR=$BASE_DIR/Plastome10x/logs

mkdir -p "$OUT_DIR" "$LOG_DIR"

# --- Main loop ---
for TAXON_DIR in "$ASSEMBLY_DIR"/*_assembly/; do
    TAXON=$(basename "$TAXON_DIR" | sed 's/_assembly//')
    FASTA=$PLASTOME_DIR/${TAXON}_cp_final.fasta
    R1=$(find "$TAXON_DIR" -maxdepth 1 -type f -name "*_R1_val_1.fq" | head -n1)
    R2=$(find "$TAXON_DIR" -maxdepth 1 -type f -name "*_R2_val_2.fq" | head -n1)

    if [[ -z "$FASTA" || -z "$R1" || -z "$R2" ]]; then
        echo "Skipping $TAXON — missing FASTA or reads"
        continue
    fi

    echo "Processing $TAXON"

    # Index plastome
    bwa index "$FASTA"

    # Map reads
    bwa mem -t 16 "$FASTA" "$R1" "$R2" | samtools view -bS - | samtools sort -@ 8 -o "$OUT_DIR/${TAXON}.bam"
    samtools index "$OUT_DIR/${TAXON}.bam"

    # Generate coverage file
    bedtools genomecov -ibam "$OUT_DIR/${TAXON}.bam" -dz > "$OUT_DIR/${TAXON}.cov"

    # Extract ≥10× regions
    awk '$3 >= 10 {print $1"\t"$2-1"\t"$2}' "$OUT_DIR/${TAXON}.cov" > "$OUT_DIR/${TAXON}_10x.bed"
    bedtools merge -i "$OUT_DIR/${TAXON}_10x.bed" > "$OUT_DIR/${TAXON}_10x_merged.bed"

    # Trim plastome to high-coverage regions
    bedtools getfasta -fi "$FASTA" -bed "$OUT_DIR/${TAXON}_10x_merged.bed" -fo "$OUT_DIR/${TAXON}_cp_10x.fasta"

    echo "Finished $TAXON"
done

# Combine all trimmed FASTAs
cat "$OUT_DIR"/*_cp_10x.fasta > "$OUT_DIR/Pritchardia_plastomes_10x_combined.fasta"

echo "ALL DONE"
date
