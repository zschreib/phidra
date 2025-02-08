import sys
from Bio import SeqIO

def remove_sequences(file1, file2, output_file):
    # Read the sequences from File 2
    seqs_to_remove = set()
    for record in SeqIO.parse(file2, "fasta"):
        seqs_to_remove.add(str(record.seq))

    # Open the output file for writing
    with open(output_file, "w") as output_handle:
        # Read sequences from File 1 and only write those not in File 2
        for record in SeqIO.parse(file1, "fasta"):
            if str(record.seq) not in seqs_to_remove:
                SeqIO.write(record, output_handle, "fasta")

    print(f"Output written to {output_file}")

if __name__ == "__main__":
    # Ensure three arguments are provided: file1, file2, output_file
    if len(sys.argv) != 4:
        print("Usage: python remove_sequences.py <file1.fa> <file2.fa> <output.fa>")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    output_file = sys.argv[3]

    remove_sequences(file1, file2, output_file)
