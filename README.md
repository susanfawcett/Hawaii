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


### Prepare paired read files and CSV to run python script SORTER2_FormatReads.py

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
export SCR="/global/scratch/users/sfawcett"

#### Conda should be activated and deactivated within the batch script to avoid problems. 
source /global/home/users/sfawcett/miniconda3/etc/profile.d/conda.sh
conda activate SORTER2

#### Include path to a single master version of scripts (these need not be modified; all necessary file names and other changes can be made to batch scripts for each project using the pipeline)
python /global/scratch/users/akostic/sorter/bin/1a.py -n Pritchardia -trim T -spades T
 
















  



       
 


