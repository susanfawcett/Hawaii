#!/bin/bash
#SBATCH --job-name=AlignLoulu
#SBATCH --account=fc_labordia
#SBATCH --partition=savio3
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/global/scratch/users/sfawcett/HybPiper/logs/AlignLoulu_%j.out
#SBATCH --error=/global/scratch/users/sfawcett/HybPiper/logs/AlignLoulu_%j.err

echo "=== Starting MAFFT alignments for Loulu nuclear loci ==="
date

# --- Environment setup ---
module load anaconda3
source /global/home/users/sfawcett/miniconda3/etc/profile.d/conda.sh
conda activate hybpiper   # MAFFT is included here; if not, replace with SORTER2

# --- Directories ---
WORKDIR=/global/scratch/users/sfawcett/HybPiper/Loulu_Data/LouluNuclearFastas
OUTDIR=$WORKDIR/AlignedFastas
mkdir -p "$OUTDIR"

cd "$WORKDIR"

# --- Alignment loop ---
echo "Aligning all .FNA files in: $WORKDIR"
for f in *.FNA; do
    base=$(basename "$f" .FNA)
    echo "  → Aligning $f ..."
    mafft --auto --thread 8 "$f" > "$OUTDIR/${base}_aln.fasta"
done

echo "=== All alignments complete ==="
date
