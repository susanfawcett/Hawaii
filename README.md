## Pipelines for the analyses of targeted enrichment data for Hawaiian Angiosperm Lineages

We include modified versions of scripts published by Jonas Mendez-Reneau, using SORTER2, and other tools.
See https://github.com/JonasMendez/SORTER2 for tutorials and detailed explanations. Scripts are modified to take advantage of multi-threading where possible, in a linux environment on UC Berkeley's Savio computing cluster. https://research-it.berkeley.edu/services-projects/high-performance-computing-savio

### List of Jonas Mendez-Reneau scripts

SORTER2_FormatReads.py
SORTER2_HaplOMiner.py
SORTER2_Processor.py
SORTER2_ProgenitorProcessor.py
SORTER2_Stage1A_TrimSPAdes.py
SORTER2_Stage1B_AssembleOrthologs.py
SORTER2_Stage2_PhaseOrthologs.py
SORTER2_Stage3_PhaseHybrids.py
example_shell.sh

### List of our Modified Scripts

1a.py
1b_pre_mafft.py 
1b_mafft.py  
1b_post_mafft.py  
2pre_bwa.py   
2bwa.py
2post_bwa.py
2end_mafft.py
3pre_mafft.py
3mafft.py
3post_mafft.py
hap1.py
hap2.py

proc_snps.py
proc1.py
proc2.py

### Install SORTER2
wget https://raw.githubusercontent.com/JonasMendez/SORTER2/refs/heads/main/SORTER2.yml

### Prepare paired read files and CSV to run python script 
Use SORTER2_FormatReads.py

#### Rename files ending in .fq with .fastq 
for f in *.fq; do
    mv "$f" "${f/.fq/.fastq}"
done

#### Rename all _1.fastq → _R1.fastq
for f in *_1.fastq; do
    mv "$f" "${f/_1.fastq/_R1.fastq}"
done

#### Rename all _2.fastq → _R2.fastq
for f in *_2.fastq; do
    mv "$f" "${f/_2.fastq/_R2.fastq}"
done

#### Remove hidden formatting from .CSV file
sed -i '1s/^\xEF\xBB\xBF//' Your_samples.csv

#### Edit .CSV file to match directory file names
sed -i -e 's/.fq/.fastq/g' -e 's/_1.fastq/_R1.fastq/g' -e 's/_2.fastq/_R2.fastq/g' Your_samples.csv

### Format Batch Script
#### See Example script. 

#### Include a path to SCR
export SCR="/path/to/directory/"

#### Conda should be activated and deactivated within the batch script to avoid problems. 
source /path/to/directory/miniconda3/etc/profile.d/conda.sh
conda activate SORTER2

### Run SORTER2 Stage 1A
Note: Include path to a single master version of python scripts. All necessary file names and other changes can be made to batch scripts for each project using the pipeline).
Batch script 
python /path/to/directory/bin/1a.py -n Projectname -trim T -spades T

### Use Hapl-O-Miner to assemble off-target reads to reference plastome
#### Reference plastomes were used for _Pritchardia pacifica_, NCBI reference NC_067842.1 and _Lysimachia clethroides_, NCBI reference  NC_064345.1, downloaded from GenBank.
Batch script prihap1.bat
Designating the trimmed paired end reads output from SORTER2 stage 1A, and the reference genome, with default threshholds -c 80 (80% reference coverage) -d 10 (10x read depth).

#### for all the folders ending in *assembly copy all files ending in *cp_final.fasta to /global/scratch/users/sfawcett/Haplominer/PritchardiaPlastomes
 
for dir in /global/scratch/users/sfawcett/SORTER2/Pritchardia/Raw/SORTER2_Pritchardia/*assembly; do
  if [ -d "$dir" ]; then
    cp "$dir"/*cp_final.fasta /global/scratch/users/sfawcett/Haplominer/PritchardiaPlastomes/ 2>/dev/null
  fi
done

#### Use filenames as headers for mafft (they otherwise will carry the name of the reference they were mapped to and are truncated)
bash plastafastomerename.sh

#### Combine all assemblies into a single fasta file 
cat *.fasta > all_loulu_plastomes_unaligned.fasta

#### Align assemblies to the reference using Mafft
Batch script MafftPritchardiaPlastome.bat
mafft --thread 40 --auto --keeplength --add "$ASSEMBLY" "$REFERENCE" > "$OUTPUT_DIR/aligned_plastomes.fasta"

#### Remove sequences with less than 50% reference coverage
python trim_by_coverage_0.5.py

#### Summary statistics can be used to evaluate different trimming strategies
python summarize_alignment.py

#### TrimAL was used to generate final alignments
For _Pritchardia_
trimal -in trimmed_alignment_50percent.fasta -out trimmed_strictplus.fasta -strictplus 
For _Lysimachia_
trimal -in trimmed_alignment_50percent.fasta -out trimmed_alignment_gappyout.fasta -gappyout 

#### Plastome Phylogenies were inferred using IQTree
Batch script IQTreeLouluPlastome.bat
  iqtree3 -s $ALIGNMENT_FILE -m TEST -bb 1000 -T AUTO -safe --prefix louluplastome

### Run SORTER2 Stage 1B
Broken into three steps to faciliate multithreading on Savio

## Extracting plastomes from WGS sequences
### Input csv file formatted into three columns with existing filename | collection ID | taxon
python SORTER2_FormatReads.py -i ArborWGSMadieae.csv

Reference plastome:
/global/scratch/users/sfawcett/WGS/Venegasia_chloroplast_genome.fasta

Reads directory:
/global/scratch/users/sfawcett/WGS/WGSMadieae

Output directory:
/global/scratch/users/sfawcett/WGS/WGS_plastomes

### Create your Conda Environment, install bioinformatics tools
conda create --name plastome_env --clone SORTER2
conda activate plastome_env
conda install -c bioconda bwa samtools bcftools fastp mafft -y

### 1. Prepare the reference**
Do this once to index the plastome:

bash
REF=/global/scratch/users/sfawcett/WGS/Venegasia_chloroplast_genome.fasta
samtools faidx "$REF"
bwa index "$REF"

### 2. Prepare a list of samples

Paired-end reads in a folder:

READS=/global/scratch/users/sfawcett/WGS/ArborWGSMadieae
OUT=/global/scratch/users/sfawcett/WGS/WGS_plastomes
mkdir -p "$OUT"/{bam,vcf,consensus,logs,trim}

ls "$READS"/*_R1.fastq | sed 's/_R1.fastq//' | xargs -n1 basename > "$OUT/samples.txt"

cat "$OUT/samples.txt"

### 3. SLURM script for mapping reads to reference, generating consensus sequences

see script map_plastomes

## 4. Combine consensus sequences after the run

bash
cat "$OUT/consensus"/*.fa > "$OUT/all_wgs_plastomes.fasta"

### 5. Align using Mafft

mafft --auto --thread 8 "$OUT/all_wgs_plastomes.fasta" > "$OUT/all_wgs_plastomes.aln.fasta"





  
  





 
















  



       
 


