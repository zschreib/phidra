import csv
import sys

#ONLY works with  (--format-mode 0) or default mmseqs output format

#Tophit by Evalue
def get_lowest_evalue_records(file_path):
    min_evalue_records = {}

    with open(file_path, mode='r', newline='') as file:
        reader = csv.reader(file, delimiter='\t')

        for row in reader:
            query_id = row[0]
            evalue = float(row[10])

            if query_id not in min_evalue_records or evalue < float(min_evalue_records[query_id]['Evalue']):
                min_evalue_records[query_id] = {
                    'Query_ID': row[0],
                    'Target_ID': row[1],
                    'Sequence_Identity': row[2],
                    'Alignment_Length': row[3],
                    'Mismatches': row[4],
                    'Gap_Openings': row[5],
                    'Query_Start': row[6],
                    'Query_End': row[7],
                    'Target_Start': row[8],
                    'Target_End': row[9],
                    'Evalue': row[10],
                    'Bit_Score': row[11]
                }

    return list(min_evalue_records.values())

#Tophit by Bitscore
def get_highest_bitscore_records(file_path):
    max_bitscore_records = {}

    with open(file_path, mode='r', newline='') as file:
        reader = csv.reader(file, delimiter='\t')

        for row in reader:
            query_id = row[0]
            bit_score = float(row[11])

            if query_id not in max_bitscore_records or bit_score > float(max_bitscore_records[query_id]['Bit_Score']):
                max_bitscore_records[query_id] = {
                    'Query_ID': row[0],
                    'Target_ID': row[1],
                    'Sequence_Identity': row[2],
                    'Alignment_Length': row[3],
                    'Mismatches': row[4],
                    'Gap_Openings': row[5],
                    'Query_Start': row[6],
                    'Query_End': row[7],
                    'Target_Start': row[8],
                    'Target_End': row[9],
                    'Evalue': row[10],
                    'Bit_Score': row[11]
                }

    return list(max_bitscore_records.values())

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python script.py <input_file> <output_evalue_file> <output_bitscore_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_evalue_file = sys.argv[2]
    output_bitscore_file = sys.argv[3]

    #Default mmseqs output headers similar to BLAST output
    headers = ['Query_ID', 'Target_ID', 'Sequence_Identity', 'Alignment_Length', 'Mismatches',
               'Gap_Openings', 'Query_Start', 'Query_End', 'Target_Start', 'Target_End', 
               'Evalue', 'Bit_Score']

    lowest_evalue_records = get_lowest_evalue_records(input_file)

    with open(output_evalue_file, mode='w', newline='') as evalue_file:
        writer = csv.DictWriter(evalue_file, fieldnames=headers, delimiter='\t')
        writer.writeheader()
        writer.writerows(lowest_evalue_records)

    highest_bitscore_records = get_highest_bitscore_records(input_file)

    with open(output_bitscore_file, mode='w', newline='') as bitscore_file:
        writer = csv.DictWriter(bitscore_file, fieldnames=headers, delimiter='\t')
        writer.writeheader()
        writer.writerows(highest_bitscore_records)

