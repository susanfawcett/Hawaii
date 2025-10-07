#!/bin/bash
#SBATCH --job-name=MaskAlign10x
#SBATCH --account=fc_labordia
#SBATCH --partition=savio3_htc
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=120G

echo "START"
date

source /global/home/users/sfawcett/miniconda3/etc/profile.d/conda.sh
conda activate SORTER2

# Activate your conda environment
conda activate SORTER2

# Directories
INPUT_DIR=/global/scratch/users/sfawcett/Haplominer/Madieae/Plastomes
TMP_DIR=$INPUT_DIR/tmp_masked
mkdir -p "$TMP_DIR"

# Loop through all assembly directories
for TAXON_DIR in "$INPUT_DIR"/*_assembly/; do
    TAXON=$(basename "$TAXON_DIR")
    FASTA=$(find "$TAXON_DIR" -maxdepth 1 -type f -name "*_cp_final.fasta" | head -n1)
    BAM=$(find "$TAXON_DIR" -maxdepth 1 -type f -name "*_cpreads_srt.bam" | head -n1)

    if [[ -z "$FASTA" || -z "$BAM" ]]; then
        echo "Skipping $TAXON_DIR: missing BAM or FASTA"
        continue
    fi

    echo "Processing $TAXON"

    # Generate BED of positions with coverage <10
    samtools depth -a "$BAM" | awk '$3 < 10 {print $1, $2-1, $2}' OFS="\t" > "$TAXON_DIR/lowcov.bed"

    # Mask low coverage positions in FASTA
    bedtools maskfasta -fi "$FASTA" -bed "$TAXON_DIR/lowcov.bed" -fo "$TMP_DIR/${TAXON}_masked.fasta"
done

echo "Combining all masked FASTAs..."
cat "$TMP_DIR"/*_masked.fasta > "$TMP_DIR/all_masked_concat.fasta"

echo "Aligning combined masked FASTAs..."
mafft --thread 40 --auto "$TMP_DIR/all_masked_concat.fasta" > "$INPUT_DIR/Madieae10x_aligned.fasta"

echo "DONE"
date
