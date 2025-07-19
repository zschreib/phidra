import argparse
import subprocess
import os
import sys
from tqdm import tqdm

# Define the bin directory
BIN_DIR = os.path.join(os.getcwd(), 'bin')
LOG_DIR = None

def set_log_dir(output_dir):
    global LOG_DIR
    LOG_DIR = output_dir
    os.makedirs(LOG_DIR, exist_ok=True)

def log_command(command):
    tqdm.write(command)
    log_file_path = os.path.join(LOG_DIR, "job_logger.txt")
    with open(log_file_path, "a") as log_file:
        log_file.write(f"{command}\n")

def run_job(command):
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        log_command(f"Error executing {command}: {e.stderr}")
        print(f"Error executing {command}: {e.stderr}", file=sys.stderr)
        sys.exit(1)

def process_input(input_fasta, subject_db, pfam_hmm_db, ida_file, output_dir, function, threads):

    log_command(f"Running {input_fasta} with the following parameters:")
    log_command(f"  Input FASTA: {input_fasta}")
    log_command(f"  Subject DB: {subject_db}")
    log_command(f"  Pfam HMM DB: {pfam_hmm_db}")
    log_command(f"  IDA File: {ida_file}")
    log_command(f"  Output Directory: {output_dir}")
    log_command(f"  Function: {function}")
    log_command(f"  Threads: {threads}")
    log_command(f"=======================================================================================")

    """Processes the input files and manages workflow steps."""
    # Define sub-directory paths
    func_dir = os.path.join(output_dir, function)
    subject_db_dir = os.path.join(func_dir, "subject_db")
    recursive_subject_db_dir = os.path.join(func_dir, "r_subject_db")

    mmseqs_results_dir = os.path.join(func_dir, "mmseqs_results", "initial_search")
    pfam_results_dir = os.path.join(func_dir, "pfam_domain_results", "initial_search")
    mmseqs_results_dir_r = os.path.join(func_dir, "mmseqs_results", "recursive_search")
    pfam_results_dir_r = os.path.join(func_dir, "pfam_domain_results", "recursive_search")

    results_dir = os.path.join(func_dir, "final_results")

    # Ensure subdirectories are created
    os.makedirs(subject_db_dir, exist_ok=True)
    os.makedirs(recursive_subject_db_dir, exist_ok=True)
    os.makedirs(mmseqs_results_dir, exist_ok=True)
    os.makedirs(pfam_results_dir, exist_ok=True)
    os.makedirs(mmseqs_results_dir_r, exist_ok=True)
    os.makedirs(pfam_results_dir_r, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    log_command(f"Running analysis for {function} sequences using {ida_file} as domain validation")
    # Define tasks with shorter path references
    tasks = [
        # First round search
        f"mmseqs createdb {subject_db} {os.path.join(subject_db_dir, f'{function}_db')}",
        f"mmseqs createindex {os.path.join(subject_db_dir, f'{function}_db')} {os.path.join(subject_db_dir, 'tmp')}",
        f"mmseqs easy-search -e 1E-3 --threads {threads} {input_fasta} {os.path.join(subject_db_dir, f'{function}_db')} {os.path.join(mmseqs_results_dir, f'{function}_results.m8')} {os.path.join(mmseqs_results_dir, 'tmp')}",
        # Top hit pull and fasta creation
        f"python {os.path.join(BIN_DIR, 'mmseqs_tophit_calculation.py')} {os.path.join(mmseqs_results_dir, f'{function}_results.m8')} {os.path.join(mmseqs_results_dir, f'{function}_TopHit_Evalue.tsv')} {os.path.join(mmseqs_results_dir, f'{function}_TopHit_Bitscore.tsv')}",
        f"python {os.path.join(BIN_DIR, 'generate_tophit_fasta.py')} {os.path.join(mmseqs_results_dir, f'{function}_TopHit_Evalue.tsv')} {input_fasta} {os.path.join(mmseqs_results_dir, f'{function}_TopHit_Evalue.fa')}",
        # Pfamscan and validation of significant domains
        f"python {os.path.join(BIN_DIR, 'pfam_scan.py')} -out {os.path.join(pfam_results_dir, f'pfam_coverage_report.tsv')} -outfmt tsv -cpu {threads} {os.path.join(mmseqs_results_dir, f'{function}_TopHit_Evalue.fa')} {pfam_hmm_db}",
        f"python {os.path.join(BIN_DIR, 'validate_pfam.py')} {os.path.join(pfam_results_dir, f'pfam_coverage_report.tsv')} {input_fasta} {ida_file} {pfam_results_dir}",
        f"python {os.path.join(BIN_DIR, 'remove_duplicates.py')} {input_fasta} {os.path.join(pfam_results_dir, f'pfam_validated_full_protein.fa')} {os.path.join(recursive_subject_db_dir, f'filtered_out_hits.fa')}",
    ]

    for i, task in enumerate(tasks, 1):
        log_command(f"[{i}/{len(tasks)}] Processing Inital Search: {task}")
        run_job(task)

    all_identified = os.path.join(recursive_subject_db_dir, f'filtered_out_hits.fa')
    if os.path.exists(all_identified) and os.path.getsize(all_identified) > 0:
        log_command(f"Not all sequences validated. Running a Recursive Search:")
        recursive_task = f"mmseqs easy-search -e 1E-3 --threads {threads} {os.path.join(recursive_subject_db_dir, f'filtered_out_hits.fa')} {os.path.join(pfam_results_dir, f'pfam_validated_full_protein.fa')} {os.path.join(mmseqs_results_dir_r, f'{function}_recursive_results.m8')} {os.path.join(mmseqs_results_dir_r, 'tmp')}"
        log_command(f"Processing Recursive Check: {recursive_task}")
        run_job(recursive_task)

    recursive_results_path = os.path.join(mmseqs_results_dir_r, f'{function}_recursive_results.m8')

    # Recusive run on results to validate possible distant hits to validated pfams/IDA matches
    if os.path.exists(recursive_results_path) and os.path.getsize(recursive_results_path) > 0:
        log_command(f"Recursive matches found. Processing {recursive_results_path}...")
        recursive_task = [

            f"python {os.path.join(BIN_DIR, 'mmseqs_tophit_calculation.py')} {recursive_results_path} {os.path.join(mmseqs_results_dir_r, f'{function}_TopHit_Evalue.tsv')} {os.path.join(mmseqs_results_dir_r, f'{function}_TopHit_Bitscore.tsv')}",
            f"python {os.path.join(BIN_DIR, 'generate_tophit_fasta.py')} {os.path.join(mmseqs_results_dir_r, f'{function}_TopHit_Evalue.tsv')} {input_fasta} {os.path.join(mmseqs_results_dir_r, f'{function}_TopHit_Evalue.fa')}",

            f"python {os.path.join(BIN_DIR, 'pfam_scan.py')} -out {os.path.join(pfam_results_dir_r, f'pfam_coverage_report.tsv')} -outfmt tsv -cpu {threads} {os.path.join(mmseqs_results_dir_r, f'{function}_TopHit_Evalue.fa')} {pfam_hmm_db}",
            f"python {os.path.join(BIN_DIR, 'validate_pfam.py')} {os.path.join(pfam_results_dir_r, f'pfam_coverage_report.tsv')} {input_fasta} {ida_file} {pfam_results_dir_r}",

            #Being merging results starting with inital search
            f"cp {os.path.join(pfam_results_dir, f'pfam_validated_full_protein.fa')} {os.path.join(pfam_results_dir, f'pfam_validated_report.tsv')} {os.path.join(pfam_results_dir, f'pfam_validated_domain_only.fa')} {results_dir}",
            #Merge recursive pfam_validated protein fasta
            f"cat {os.path.join(pfam_results_dir_r, f'pfam_validated_full_protein.fa')} >> {os.path.join(results_dir, f'pfam_validated_full_protein.fa')}",
            #Merge recursive pfam_domain fasta
            f"cat {os.path.join(pfam_results_dir_r, f'pfam_validated_domain_only.fa')} >> {os.path.join(results_dir, f'pfam_validated_domain_only.fa')}",
            #Merge pfam validation report
            f"tail -n +2 {os.path.join(pfam_results_dir_r, f'pfam_validated_report.tsv')} >> {os.path.join(results_dir, f'pfam_validated_report.tsv')}"

       ]
    else:
        log_command(f"Warning: {recursive_results_path} is empty, no matches found. Skipping recursive pfam scan and IDA validation.")
        recursive_task = [
            f"printf 'No recursive matches found.' > {recursive_results_path}",
            f"printf 'No recursive matches found.' > {os.path.join(pfam_results_dir_r, f'pfam_coverage_report.tsv')}",
            f"mv {os.path.join(pfam_results_dir, f'pfam_validated_full_protein.fa')} {os.path.join(pfam_results_dir, f'pfam_validated_report.tsv')} {os.path.join(pfam_results_dir, f'pfam_validated_domain_only.fa')} {results_dir}"
        ]

    for i, task in enumerate(recursive_task, 1):
        log_command(f"[{i}/{len(recursive_task)}] Processing Recursive Run: {task}")
        run_job(task)

    # Clean up tasks
    clean_up = [
        f"rm -rf {os.path.join(mmseqs_results_dir, 'tmp')}",
        f"rm -rf {os.path.join(mmseqs_results_dir_r, 'tmp')}",
        f"rm -rf {subject_db_dir}",
        f"rm -rf {recursive_subject_db_dir}",
    ]

    for i, task in enumerate(clean_up, 1):
        log_command(f"[{i}/{len(clean_up)}] Processing Clean up: {task}")
        run_job(task)

    log_command(f"Job has successfully completed!")

def main():
    parser = argparse.ArgumentParser(
        description="Identifies homologous proteins and associated Pfam domains from input protein sequences, while comparing against InterPro Domain Architectures to analyze domain-level similarities and functional relationships.",
        add_help=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Help arguments
    help_args = parser.add_argument_group("Help")
    help_args.add_argument("--help", "-h", action="help", help="Description of tool and usage")

    # Required arguments
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("--input_fasta", "-i", required=True, help="Input protein (aa) FASTA file")
    required_args.add_argument("--subject_db", "-db", required=True, help="Protein subject database FASTA file")
    required_args.add_argument("--pfam_hmm_db", "-pfam", required=True, help="Pfam subject mmseqs DB file")
    required_args.add_argument("--ida_file", "-ida", required=True, help="IDA file for processing. Must be in Interpro IDA TSV format")
    required_args.add_argument("--function", "-f", required=True, help="Function of the protein of interest")
    required_args.add_argument("--output_dir", "-o", required=True, help="Directory for output results of that project")
    
    # Optional arguments group
    optional_args = parser.add_argument_group("Optional Arguments")
    optional_args.add_argument("--threads", "-t", type=int, default=2, help="Number of threads to use")
    optional_args.add_argument("--version", "-v", action="version", version="PHIDRA v1.0", help="Show program's version number and exit")

    args = parser.parse_args()

    protein_analysis = os.path.join(args.output_dir, args.function)

    # Check if output directory exists
    if os.path.exists(protein_analysis):
        print(f"Output project for protein '{protein_analysis}' already exists. To avoid overriding data, please rename or delete it and try again.")
        sys.exit(1)
    else:
        os.makedirs(protein_analysis)
        set_log_dir(protein_analysis)

    # Process inputs
    process_input(args.input_fasta, args.subject_db, args.pfam_hmm_db, args.ida_file, args.output_dir, args.function, args.threads)

if __name__ == "__main__":
    main()
