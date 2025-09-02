#!/bin/bash

# Plastome mapping → consensus pipeline
# Adapted for Savio, unzipped reads, single/paired

#SBATCH --job-name=wgs_plastomes
#SBATCH --account=fc_labordia
#SBATCH --partition=savio3
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-93  

echo
echo "START"; date
echo

# Load modules and activate conda environment
module load anaconda3 bwa samtools bcftools fastp
source activate plastome_env

# Paths
REF=/global/scratch/users/sfawcett/WGS/Venegasia_chloroplast_genome.fasta
READS=/global/scratch/users/sfawcett/WGS/ArborWGSMadieae
OUT=/global/scratch/users/sfawcett/WGS/WGS_plastomes

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$OUT/samples.txt")
LOG="$OUT/logs/${SAMPLE}.log"

echo "[$(date)] Starting $SAMPLE" > "$LOG"

# Assign input reads
R1="$READS/${SAMPLE}_R1.fastq"
R2="$READS/${SAMPLE}_R2.fastq"

# Check if R2 exists (single-end vs paired-end)
if [[ ! -f "$R2" ]]; then
    R2=""
fi

# Trim reads
if [[ -n "$R2" ]]; then
    fastp -i "$R1" -I "$R2" \
          -o "$OUT/trim/${SAMPLE}_R1.trim.fq" \
          -O "$OUT/trim/${SAMPLE}_R2.trim.fq" \
          -w $SLURM_CPUS_PER_TASK \
          -h "$OUT/logs/${SAMPLE}_fastp.html" \
          -j "$OUT/logs/${SAMPLE}_fastp.json" >> "$LOG" 2>&1
    R1="$OUT/trim/${SAMPLE}_R1.trim.fq"
    R2="$OUT/trim/${SAMPLE}_R2.trim.fq"
else
    fastp -i "$R1" \
          -o "$OUT/trim/${SAMPLE}_R1.trim.fq" \
          -w $SLURM_CPUS_PER_TASK \
          -h "$OUT/logs/${SAMPLE}_fastp.html" \
          -j "$OUT/logs/${SAMPLE}_fastp.json" >> "$LOG" 2>&1
    R1="$OUT/trim/${SAMPLE}_R1.trim.fq"
fi

# Map reads
if [[ -n "$R2" ]]; then
    bwa mem -t $SLURM_CPUS_PER_TASK "$REF" "$R1" "$R2" \
        | samtools view -b - \
        | samtools sort -@ $SLURM_CPUS_PER_TASK -o "$OUT/bam/${SAMPLE}.sorted.bam" - >> "$LOG" 2>&1
else
    bwa mem -t $SLURM_CPUS_PER_TASK "$REF" "$R1" \
        | samtools view -b - \
        | samtools sort -@ $SLURM_CPUS_PER_TASK -o "$OUT/bam/${SAMPLE}.sorted.bam" - >> "$LOG" 2>&1
fi

samtools index "$OUT/bam/${SAMPLE}.sorted.bam"

# Filter for quality
samtools view -@ $SLURM_CPUS_PER_TASK -b -q 30 -f 0x2 "$OUT/bam/${SAMPLE}.sorted.bam" > "$OUT/bam/${SAMPLE}.filtered.bam"
samtools index "$OUT/bam/${SAMPLE}.filtered.bam"

# Call variants
bcftools mpileup -f "$REF" "$OUT/bam/${SAMPLE}.filtered.bam" -Ou \
    | bcftools call -mv -Oz -o "$OUT/vcf/${SAMPLE}.vcf.gz"
bcftools index "$OUT/vcf/${SAMPLE}.vcf.gz"

# Mask low coverage sites (<5x)
samtools depth -a "$OUT/bam/${SAMPLE}.filtered.bam" > "$OUT/bam/${SAMPLE}.depth.txt"
awk '$3<5{print $1"\t"$2-1"\t"$2}' "$OUT/bam/${SAMPLE}.depth.txt" > "$OUT/bam/${SAMPLE}.lowcov.bed"

# Build consensus sequence
bcftools consensus -f "$REF" -m "$OUT/bam/${SAMPLE}.lowcov.bed" "$OUT/vcf/${SAMPLE}.vcf.gz" > "$OUT/consensus/${SAMPLE}.fa"

echo "[$(date)] Finished $SAMPLE" >> "$LOG"
echo "DONE"; date

