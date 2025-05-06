from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

alignment = AlignIO.read("LouluPlastomes5May2025.fasta", "fasta")
min_coverage = 0.5  # 50% threshold
alignment_length = alignment.get_alignment_length()

filtered = []
for record in alignment:
    ungapped = record.seq.count("-") + record.seq.count("N")
    coverage = (alignment_length - ungapped) / alignment_length
    if coverage >= min_coverage:
        filtered.append(record)

print(f"Kept {len(filtered)} of {len(alignment)} sequences.")

output_file = "trimmed_alignment_50percent.fasta"
AlignIO.write(MultipleSeqAlignment(filtered), output_file, "fasta")
print(f"Trimmed alignment saved to {output_file}")
