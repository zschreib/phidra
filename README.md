#  PHIDRA ![Version](https://img.shields.io/badge/version-1.0_stable-green)
**P**rotein **H**omology **I**dentification via **D**omain-**R**elated **A**rchitecture

A simple way to search and validate identifed Pfam domains of interest against a curated InterProScan Domain Architecture (IDA) file to check whether or not your proteins match a domain composition found in the Interpro Database.

## Description

A Python-based package of scripts that performs an initial homology search using `MMseqs2` to identify targeted proteins of interest in small or large datasets. Top hits are searched against the Pfam database using `pfam_scan`, and verified domains are checked and compared against a custom input InterProScan Domain Architecture (IDA) file relative to your target protein of interest. 
A recursive search is then performed using the full length proteins with validated IDAs as the subject database and the original input as the query, with initial matches filtered out. This process captures potentially more distant proteins that may have been missed in the initial homology search but are functionally relevant.

## Getting Started

### Dependencies

* `phidra` can be ran either on a machine with a properly set up python environment or by creating a custom [Conda](https://docs.conda.io/en/latest/) environment if you do not wish to change your current setup (recommended).
* [MMseqs2](https://github.com/soedinglab/MMseqs2) required for inital and recursive homology search. 
* [HMMER](http://hmmer.org/documentation.html) required for pfam database creation and protein domain identifcation through `pfam_scan`.
* [Python](https://www.python.org/downloads/) >= 3.8 
### Installation

* Conda setup (recommended) [Help](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
```
#Create environment with >= Python 3.8
conda create --name phidra python=3.8

#Activate environment
conda activate phidra

#Install hmmer/mmseqs2 (Required)
conda install -c conda-forge -c bioconda mmseqs2
conda install bioconda::hmmer

#Clone tool into working directory
git clone https://github.com/zschreib/phidra

cd phidra

#Grabs required python packages
pip install -r requirements.txt

```
* Python setup >= Python 3.8, requires MMseqs2 and HMMER to be already set up.
```
git clone https://github.com/zschreib/phidra
cd phidra
pip install -r requirements.txt
```

* Tool Check. If you recieve any errors here do not proceede until fixed.
```
mmseqs -h
hmmscan -h
python phidra_run.py -v
```
### Setting up Pfam Database (Optional)
* If already have Pfam-hmm formatted database you can skip this step. 

```
mkdir pfam_database
cd pfam_database

wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

gunzip -c Pfam-A.hmm.dat.gz > Pfam-A.hmm.dat
gunzip -c Pfam-A.hmm.gz > Pfam-A.hmm
rm Pfam-A.hmm.gz Pfam-A.hmm.dat.gz

hmmpress Pfam-A.hmm
```

## Creating an InterProScan Domain Architecture (IDA) file

1. Identify essential domains
* For DNA polymerase I, we have identified PF00476 (DNA_pol_A) as core domain 
* Search domain on InterPro: `https://www.ebi.ac.uk/interpro/entry/pfam/PF00476/domain_architecture/`

2. Download domain architecture data in TSV format for all or selected domain combinations\
File should include:
* IDA ID (unique hash)
* Domain combinations
* Protein counts within IDA
* Representative sequence
* Representative length
* Domain positions/coordinates

3. Sample IDA file in TSV format:
```
IDA ID	IDA Text	Unique Proteins	Representative Accession	Representative Length	Representative Domains
82911c7e8cf0ed5121595d5944b2a2f9a2c4f49e	PF02739:IPR020046-PF01367:IPR020045-PF01612:IPR002562-PF00476:IPR001098	22766	P00582	928	PF02739{5_3_exonuc_N}:IPR020046{5-3_exonucl_a-hlix_arch_N}[9-170],PF01367{5_3_exonuc}:IPR020045{DNA_polI_H3TH}[171-272],PF01612{DNA_pol_A_exo1}:IPR002562{3'-5'_exonuclease_dom}[330-516],PF00476{DNA_pol_A}:IPR001098{DNA-dir_DNA_pol_A_palm_dom}[551-925]
d4033548d92ed83469d80faa39cea59a849060a8	PF00476:IPR001098	18435	P00581	704	PF00476{DNA_pol_A}:IPR001098{DNA-dir_DNA_pol_A_palm_dom}[333-701]
16983c3a921f1409678537d9977eca4ed176e7c2	PF02739:IPR020046-PF01367:IPR020045-PF22619:IPR054690-PF00476:IPR001098	10607	Q04957	877	PF02739{5_3_exonuc_N}:IPR020046{5-3_exonucl_a-hlix_arch_N}[4-170],PF01367{5_3_exonuc}:IPR020045{DNA_polI_H3TH}[172-266],PF22619{DNA_polI_exo1}:IPR054690{DNA_polI_exonuclease}[318-457],PF00476{DNA_pol_A}:IPR001098{DNA-dir_DNA_pol_A_palm_dom}[498-875]
ce54c2b3dc9b223aa1f4992353c0972accb51c74	PF02739:IPR020046-PF01367:IPR020045-PF00476:IPR001098	6038	O84500	866	PF02739{5_3_exonuc_N}:IPR020046{5-3_exonucl_a-hlix_arch_N}[3-166],PF01367{5_3_exonuc}:IPR020045{DNA_polI_H3TH}[167-259],PF00476{DNA_pol_A}:IPR001098{DNA-dir_DNA_pol_A_palm_dom}[495-865]
f54dc0def957e931720f1460942192debd6092ca	PF01612:IPR002562-PF00476:IPR001098	4576	Q05254	595	PF01612{DNA_pol_A_exo1}:IPR002562{3'-5'_exonuclease_dom}[19-210],PF00476{DNA_pol_A}:IPR001098{DNA-dir_DNA_pol_A_palm_dom}[243-587]
```
4. Analyze architectures


Example: Core domain (PF00476) appears in different combinations

Common associated domains:
* PF02739 (5' exonuclease N-terminal)
* PF01367 (5' exonuclease)
* PF01612 (DNA polymerase A exonuclease)

5. Use for validation

* Compare unknown or known sequences against known architectures
* Verify domain organization matches canonical patterns
* Check domain order and spacing is biologically plausible (if desired)

Domain architectures help validate protein predictions by ensuring essential functional domains are present in correct arrangements.


## Running the tool
```
usage: phidra_run.py [--help] --input_fasta INPUT_FASTA --subject_db SUBJECT_DB --pfam_hmm_db PFAM_HMM_DB
--ida_file IDA_FILE --function FUNCTION --output_dir OUTPUT_DIR [--threads THREADS] [--version]

Identifies homologous proteins and associated Pfam domains from input protein sequences,
while comparing against InterPro Domain Architectures to analyze domain-level similarities and functional relationships.

Help:
  --help, -h            Description of tool and usage

Required arguments:
  --input_fasta INPUT_FASTA, -i INPUT_FASTA
                        Input protein (aa) FASTA file (default: None)
  --subject_db SUBJECT_DB, -db SUBJECT_DB
                        Protein subject database FASTA file (default: None)
  --pfam_hmm_db PFAM_HMM_DB, -pfam PFAM_HMM_DB
                        Pfam subject mmseqs DB file (default: None)
  --ida_file IDA_FILE, -ida IDA_FILE
                        IDA file for processing. Must be in Interpro IDA TSV format (default: None)
  --function FUNCTION, -f FUNCTION
                        Function of the protein of interest (default: None)
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        Directory for output results of that project (default: None)

Optional Arguments:
  --threads THREADS, -t THREADS
                        Number of threads to use (default: 2)
  --version, -v         Show program's version number and exit
```


* Example run using the test data to search for DNA polymerase A
```
python phidra_run.py -i test_data/test_query/4_seq_test.fa -db test_data/test_subject/pola_16_k12_ref.fa -pfam pfam_database/pfam_hmm -ida test_data/ida_files/polA_IDA_list.tsv -f POLA -o test_results
```

## Results 

```
test_results/
└── POLA
    ├── final_results
    │   ├── pfam_validated_domain_only.fa
    │   ├── pfam_validated_full_protein.fa
    │   └── pfam_validated_report.tsv
    ├── mmseqs_results
    │   ├── initial_search
    │   │   ├── POLA_results.m8
    │   │   ├── POLA_TopHit_Bitscore.tsv
    │   │   ├── POLA_TopHit_Evalue.fa
    │   │   └── POLA_TopHit_Evalue.tsv
    │   └── recursive_search
    │       └── POLA_recursive_results.m8
    └── pfam_domain_results
        ├── initial_search
        │   └── pfam_coverage_report.tsv
        └── recursive_search
            └── pfam_coverage_report.tsv
```

ID: ENA_MF036692_MF036692_1_29921_27213_46 was apart of the DNA polymerase family B, did not have valid domains and therefore was not included in the final pfam_validated_reports

## Authors

Contributors names and contact info

zschreib@udel.edu

## License

This project is licensed under the GNU General Public License v3.0 see the LICENSE file for details

## Acknowledgments

Inspiration, code snippets, etc.
* [pfam_scan](https://github.com/aziele/pfam_scan)
