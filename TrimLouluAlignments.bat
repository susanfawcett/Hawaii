#!/bin/bash
#SBATCH --job-name=TrimLouluSeq
#SBATCH --account=fc_labordia
#SBATCH --partition=savio3
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

module load anaconda3
source /global/home/users/sfawcett/miniconda3/etc/profile.d/conda.sh
conda activate hybpiper

WORKDIR=/global/scratch/users/sfawcett/HybPiper/Loulu_Data/LouluSuperContigs/AlignedSuperContigs
OUTDIR=/global/scratch/users/sfawcett/HybPiper/Loulu_Data/LouluSuperContigs/TrimmedAlignedSuperContigs
mkdir -p "$OUTDIR"

TRIMAL=/global/home/users/sfawcett/miniconda3/envs/hybpiper/bin/trimal
LOGFILE="$OUTDIR/trim_report.tsv"

# Header for report
echo -e "Locus\tOriginalLength\tTrimmedLength\tRemovedSequences\tTrimAlCriteria" > "$LOGFILE"

cd "$WORKDIR"

for f in *_aln.fasta; do
    base=$(basename "$f" _aln.fasta)
    tmp_upper="$OUTDIR/${base}_UPPER.tmp"
    out="$OUTDIR/${base}_trim.fasta"

    # Uppercase sequences
    awk '{
        if ($0 ~ /^>/) print $0;
        else print toupper($0);
    }' "$f" > "$tmp_upper"

    # Original alignment length (number of columns)
    ORIG_LEN=$(grep -v '^>' "$tmp_upper" | awk '{print length}' | head -n 1)

    # Run trimAl, capture warnings
    WARNINGS=$($TRIMAL -in "$tmp_upper" -out "$out" -automated1 2>&1 | grep "composed only by gaps" | sed 's/WARNING: Removing sequence //')

    # Trimmed alignment length
    TRIM_LEN=$(grep -v '^>' "$out" | awk '{print length}' | head -n 1)

    # Record locus, lengths, removed sequences, criteria
    if [ -n "$WARNINGS" ]; then
        REMOVED="$WARNINGS"
    else
        REMOVED="None"
    fi

    echo -e "${base}\t${ORIG_LEN}\t${TRIM_LEN}\t${REMOVED}\t-automated1" >> "$LOGFILE"

    rm "$tmp_upper"
done

echo "Done. Full trimming report saved in $LOGFILE"

