## Pipelines for the analyses of targeted enrichment data for Hawaiian Angiosperm Lineages

### HybPiper used for assemblies
Following tutorials on https://github.com/mossmatters/HybPiper


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
    SORTER2_FormatReads.py

    example: Pritchardia_waialealeana___Kauai_Horn_4950_ILM166_R1.fastq,Horn4950Kauai,Pri.waialealeana

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
Note: Include path to a single master version of python scripts. All necessary file names and other changes can be made to batch scripts for each project using the pipeline.
Batch script 
    python /path/to/directory/bin/1a.py -n Projectname -trim T -spades T

### Run SORTER2 Stage 1B




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

For _Geniostoma_ 99.99% average coverage, one sequence removed for failure Gen.tinifoliumtin_OppenheimerH62005Maui

    

#### Plastome Phylogenies were inferred using IQTree
Batch script IQTreeLouluPlastome.bat

      iqtree3 -s $ALIGNMENT_FILE -m TEST -bb 1000 -T AUTO -safe --prefix louluplastome

### Run SORTER2 Stage 1B
Broken into three steps to faciliate multithreading on Savio



## Extracting plastomes from WGS sequences using GetOrganelle
### Input csv file formatted into three columns with existing filename | collection ID | taxon
    python SORTER2_FormatReads.py -i ArborWGSMadieae.csv

Reference plastome:
/global/scratch/users/sfawcett/WGS/Venegasia_chloroplast_genome.fasta

Reads directory:
/global/scratch/users/sfawcett/WGS/WGSMadieae

Output directory:
/global/scratch/users/sfawcett/WGS/WGS_plastomes

### Create your Conda Environment, install bioinformatics tools
Follow steps on https://github.com/Kinggerm/GetOrganelle

### GetOrganelle 
See shell script getorganelle.sh

Jian-Jun Jin*, Wen-Bin Yu*, Jun-Bo Yang, Yu Song, Claude W. dePamphilis, Ting-Shuang Yi, De-Zhu Li. GetOrganelle: a fast and versatile toolkit for accurate de novo assembly of organelle genomes. Genome Biology 21, 241 (2020). https://doi.org/10.1186/s13059-020-02154-5

Venegesia carpesioides was used as a seed, default settings for embryophyta were used.

#### Use best Assemblies for single folder of FASTA files
    #!/usr/bin/env bash
    
    #### run from: /global/scratch/users/sfawcett/GetOrganelle/GetOrganelleOutput/

    OUTDIR=/global/scratch/users/sfawcett/GetOrganelle/GetOrganelleAssemblies
    mkdir -p "$OUTDIR"
    MAPFILE="$OUTDIR/fasta_name_mapping.txt"
    echo -e "old_name\tnew_name" > "$MAPFILE"

    for d in */ ; do
    sample=$(basename "$d")
    log="$d/get_org.log.txt"
    [ -f "$log" ] || { echo "No log for $sample, skipping"; continue; }

    #### find the line that says "Writing GRAPH to ...selected_graph.gfa"
    graph_line=$(grep "Writing GRAPH to" "$log")
    #### extract path number from the line before it (the last PATH written)
    #### assuming PATH lines appear right before GRAPH line
    path_line_num=$(grep -n "Writing PATH" "$log" | tail -n1 | cut -d: -f1)
    path_file=$(sed -n "${path_line_num}p" "$log" | awk '{print $NF}')

    if [ ! -f "$path_file" ]; then
        echo "Selected path not found for $sample, skipping"
        continue
    fi

    #### copy to new directory with folder name
    new_fasta="$OUTDIR/${sample}.fasta"
    cp "$path_file" "$new_fasta"

    #### record mapping
    echo -e "$(basename "$path_file")\t${sample}.fasta" >> "$MAPFILE"
    done

#### Errors or poor sequences resulted in no assemblies for:
    BGB668_DUBAUTIApauciflorula
    PWBH925_DUBAUTIAlinearisLIN

#### Samples removed for non-circular plastomes, low sequence length, etc.

    BGB667_DUBAUTIAimbricataIMB.fasta  	26	1061   
    BGB603_JENSIAyosemitana.fasta      	2 	108195 
    GetOrganellePlastomes.fasta        	18	151864 
    BGB675_DUBAUTIAlatifolia.fasta     	3 	18369  
    BGB509b_KYHOSIAbolanderi.fasta     	2 	18390  
    BGB664_DUBAUTIAreticulata.fasta    	2 	18393  
    BGB538_DEINANDRApentactis.fasta    	1 	189069 
    BGB670_DUBAUTIAraillardioides.fasta	22	2506   
    KNFGsn1982_DUBAUTIAlinearisLIN.fast	12	3226   
    BGB662_DUBAUTIAlaxaLAX.fasta       	34	670    
    BGB1231_DEINANDRApaniculataCRU.fast	4 	9476   

#### Rewrite fasta headers to match file names from SubsetGetOrganelleAssemblies directory
    for f in *.fasta; do base=$(basename "$f" .fasta); awk -v name="$base" '/^>/{print ">"name; next} {print}' "$f" > tmp && mv tmp "$f"; done

#### Concatenate for Mafft
    cat *.fasta > /global/scratch/users/sfawcett/GetOrganelle/SubsetGetOrganelleAssemblies/SubsetGetOrganellePlastomes.fasta

#### Assess sequences manually in Geneious 

Erroneous base calling in Haplominer introduced directional bias resulting in artefactual similarity reflecting method rather than phylogenetic relationships. The BAM and FASTA files from the Haplominer assembly were used to filter and mask plastomes at 10x coverage for further comparison with GetOrganelle results, and aligned using Mafft.

See MaskAlign10x.bat

These were then aligned with the GetOrganelle plastomes, which were more complete and had higher depth of coverage.



#### 



## Analyses of Nuclear data for Madieae
#### Combine HybPhaser exon only phased AMA_ASC phylogeny with diploids

    mkdir -p MergedSwordsDiploids

    for f in out_align_exons_192/*.fasta; do
        fname=$(basename "$f")
        cat "$f" Silversword_AMA_ASC/"$fname" 2>/dev/null > MergedSwordsDiploids/"$fname"
    done


    
    





  
  





 
















  



       
 


