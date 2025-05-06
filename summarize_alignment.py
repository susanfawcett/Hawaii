from Bio import AlignIO

#alignment = AlignIO.read("LouluPlastomes5May2025.fasta", "fasta")
#alignment = AlignIO.read("trimmed_alignment_gapsremoved0.5.fasta", "fasta")
#alignment = AlignIO.read("trimmed_alignment_nogaps.fasta", "fasta")
#alignment = AlignIO.read("aggro_trimmed.fasta", "fasta")
#alignment = AlignIO.read("trimmed_gt0.9.fasta", "fasta")
alignment = AlignIO.read("trimmed_strictplus.fasta", "fasta")

n_sequences = len(alignment)
alignment_length = alignment.get_alignment_length()

print(f"Number of sequences: {n_sequences}")
print(f"Alignment length: {alignment_length}")

# Count % of non-gap characters per sequence
for record in alignment:
    ungapped = record.seq.count("-") + record.seq.count("N")
    coverage = 100 * (alignment_length - ungapped) / alignment_length
    print(f"{record.id}\t{coverage:.2f}% coverage")
