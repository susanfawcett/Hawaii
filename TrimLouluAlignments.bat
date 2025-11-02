#!/bin/bash
#SBATCH --job-name=TrimLouluSeq
#SBATCH --account=fc_labordia
#SBATCH --partition=savio3
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G


echo "=== SLURM job starting ==="
date

module load anaconda3
source /global/home/users/sfawcett/miniconda3/etc/profile.d/conda.sh
conda activate hybpiper

echo "Conda env: $CONDA_DEFAULT_ENV"
echo "trimal is at:" $(which trimal)
trimal -h | head -n 3

WORKDIR=/global/scratch/users/sfawcett/HybPiper/Loulu_Data/LouluNuclearFastas/AlignedFastas
OUTDIR=/global/scratch/users/sfawcett/HybPiper/Loulu_Data/LouluNuclearFastas/TrimmedAlignedFastas
mkdir -p "$OUTDIR"

echo "Files in WORKDIR:"
ls -l "$WORKDIR"

# --- Environment ---
module load anaconda3
source /global/home/users/sfawcett/miniconda3/etc/profile.d/conda.sh
conda activate hybpiper

TRIMAL=/global/home/users/sfawcett/miniconda3/envs/hybpiper/bin/trimal
echo "Using trimAl: $TRIMAL"
$TRIMAL -h | head -n 3

cd "$WORKDIR"

for f in *_aln.fasta; do
    base=$(basename "$f" _aln.fasta)
    echo "Trimming $f ..."
    $TRIMAL -in "$f" -out "$OUTDIR/${base}_trim.fasta" -automated1 || { echo "ERROR: trimal failed on $f"; exit 1; }
done

echo "=== All trimAl operations complete ==="
date
