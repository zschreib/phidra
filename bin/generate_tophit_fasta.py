import sys
from Bio import SeqIO

def generate_tophit_fasta(tophit_file, fasta_file, output_fasta):

    query_ids = set()
    with open(tophit_file, "r") as file:
        next(file)
        for line in file:
            fields = line.strip().split("\t")
            query_id = fields[0]
            query_ids.add(query_id)

    print(query_ids)

    # Parse the original query FASTA and create a dictionary of sequences
    sequences = {record.id: record for record in SeqIO.parse(fasta_file, "fasta")}
    print(sequences)

    # Write out the matching sequences to a new FASTA file
    with open(output_fasta, "w") as out_fasta:
        for seq_id in query_ids:
            if seq_id in sequences:
                SeqIO.write(sequences[seq_id], out_fasta, "fasta")

    print(f"Created {output_fasta}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python generate_tophit_fasta.py <TopHit_Evalue.tsv> <Original_FASTA_File> <Output_FASTA_File>")
        sys.exit(1)

    tophit_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_fasta = sys.argv[3]

    generate_tophit_fasta(tophit_file, fasta_file, output_fasta)
