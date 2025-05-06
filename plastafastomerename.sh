# Move into the folder with your original FASTA files
cd /global/scratch/users/sfawcett/Haplominer/PritchardiaPlastomes

# Make a new directory for the renamed files
mkdir -p ../PritchardiaPlastomesRenamed

# Loop over each FASTA file
for file in *.fasta; do
    basename_no_ext="${file%.fasta}"  # strip .fasta from filename

    # Write a new file into the new directory
    {
        echo ">$basename_no_ext"      # new header
        tail -n +2 "$file"             # rest of the file (sequence only)
    } > "../PritchardiaPlastomesRenamed/$file"
done
