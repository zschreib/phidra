import sys
import csv
from collections import defaultdict, Counter
from Bio import SeqIO

def process_pfam_scan(pfam_scan_file):
    aggregated_data = defaultdict(lambda: {
        "Pfam_IDs": [],
        "Pfam_Functions": [],
        "Env_Starts": [],
        "Env_Ends": [],
        "IDA_ID": None,
        "IDA_Text": None,
        "Representative_Domains": None
    })

    with open(pfam_scan_file, 'r') as pfam_file:
        reader = csv.DictReader(pfam_file, delimiter='\t')

        for row in reader:
            if row['significance'] == '1':
                seq_id = row['seq_id']
                pfam_id = row['hmm_acc'].split('.')[0]
                pfam_function = row['hmm_name']
                env_start = int(row['env_start'])
                env_end = int(row['env_end'])

                # Aggregate data for the seq_id
                aggregated_data[seq_id]["Pfam_IDs"].append(pfam_id)
                aggregated_data[seq_id]["Pfam_Functions"].append(pfam_function)
                aggregated_data[seq_id]["Env_Starts"].append(env_start)
                aggregated_data[seq_id]["Env_Ends"].append(env_end)

    return aggregated_data

def process_ida_text(ida_text):
    """Converts IDA Text to a standardized Pfam chain format by removing IPR IDs."""
    pfam_list = []
    for entry in ida_text.split('-'):
        # Extract only the Pfam part before any ':' character
        pfam = entry.split(':')[0]
        pfam_list.append(pfam)
    # Join with '|' to match the format of the Pfam_ID chain
    return "|".join(pfam_list)

def compare_chains(chain1, chain2):
    """Compares two Pfam chains to check if they have the same Pfam occurrences."""
    return Counter(chain1.split('|')) == Counter(chain2.split('|'))

def filter_with_ida_file(aggregated_data, ida_file):
    """
    Filters the aggregated data based on Pfam chain matches with the IDA file.
    Adds IDA metadata for matched entries.
    """
    valid_pfam_chains = {}

    # Read and process the IDA file
    with open(ida_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            raw_ida_text = row['IDA Text']  # Raw Pfam chain from the IDA file
            processed_chain = process_ida_text(raw_ida_text)  # Process into a standardized chain
            valid_pfam_chains[processed_chain] = {
                "IDA_ID": row['IDA ID'],  # Assume 'IDA_ID' column exists
                "IDA_Text": raw_ida_text,
                "Representative_Domains": row['Representative Domains']
            }

    # Filter aggregated data
    filtered_data = {}
    for seq_id, data in aggregated_data.items():
        pfam_chain = "|".join(data["Pfam_IDs"])  # Construct the Pfam chain
        for valid_chain, ida_metadata in valid_pfam_chains.items():
            if compare_chains(pfam_chain, valid_chain):
                filtered_data[seq_id] = data
                filtered_data[seq_id].update(ida_metadata)
                break

    print(f"Filtered data size: {len(filtered_data)} sequences remain after matching IDA chains.")
    return filtered_data

def write_aggregated_report(aggregated_data, output_report):
    with open(output_report, 'w') as output_handle:
        writer = csv.DictWriter(
            output_handle,
            fieldnames=["Query_ID", "Pfam_IDs", "Env_Starts", "Env_Ends", "Pfam_Functions", "IDA_ID", "IDA_Text", "Representative_Domains"],
            delimiter='\t'
        )
        writer.writeheader()
        for seq_id, data in aggregated_data.items():
            writer.writerow({
                "Query_ID": seq_id,
                "Pfam_IDs": "|".join(data["Pfam_IDs"]),
                "Env_Starts": "|".join(map(str, data["Env_Starts"])),
                "Env_Ends": "|".join(map(str, data["Env_Ends"])),
                "Pfam_Functions": "|".join(data["Pfam_Functions"]),
                "IDA_ID": data.get("IDA_ID", ""),
                "IDA_Text": data.get("IDA_Text", ""),
                "Representative_Domains": data.get("Representative_Domains", "")
            })
    print(f"Aggregated report saved to: {output_report}")

def extract_sequences(fasta_file, aggregated_data, output_fasta):
    with open(output_fasta, 'w') as output_handle:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_id = record.id

            if seq_id in aggregated_data:
                data = aggregated_data[seq_id]

                for pfam_function, env_start, env_end in zip(data["Pfam_Functions"], data["Env_Starts"], data["Env_Ends"]):
                    # Extract the relevant region (env_start and env_end are 1-based inclusive)
                    extracted_seq = record.seq[env_start - 1:env_end]

                    new_header = f"{record.id}_{pfam_function}"
                    output_handle.write(f">{new_header}\n{extracted_seq}\n")

    print(f"Extracted sequences saved to: {output_fasta}")

def output_full_proteins(fasta_file, aggregated_data, output_fasta):
    with open(output_fasta, 'w') as output_handle:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_id = record.id

            if seq_id in aggregated_data:
                output_handle.write(f">{record.id}\n{record.seq}\n")

    print(f"Full protein sequences saved to: {output_fasta}")

def main(pfam_scan_file, fasta_file, ida_file, output_location):

    # Ensure the output location ends with a slash
    if not output_location.endswith('/'):
        output_location += '/'

    # Define output file paths
    aggregated_report = f"{output_location}pfam_validated_report.tsv"
    extracted_sequences = f"{output_location}pfam_validated_domain_only.fa"
    full_protein_sequences = f"{output_location}pfam_validated_full_protein.fa"

    aggregated_data = process_pfam_scan(pfam_scan_file)
    filtered_data = filter_with_ida_file(aggregated_data, ida_file)
    
    write_aggregated_report(filtered_data, output_report=aggregated_report)
    extract_sequences(fasta_file, filtered_data, output_fasta=extracted_sequences)
    output_full_proteins(fasta_file, filtered_data, output_fasta=full_protein_sequences)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <pfam_scan_output.tsv> <protein_fasta_file> <ida_file.tsv> <output_location_dir>")
        sys.exit(1)

    pfam_scan_file = sys.argv[1]
    fasta_file = sys.argv[2]
    ida_file = sys.argv[3]
    output_location = sys.argv[4]

    main(pfam_scan_file, fasta_file, ida_file, output_location)
