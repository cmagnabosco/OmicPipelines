# OmicPipelines
A collection of codes to run meta-omic analyses

## Requirements
Python 2.7, Biopython, Networkx

PRICE assembler (http://derisilab.ucsf.edu/software/price/) set to PATH

BLAST (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) set to PATH 

PhyML (https://code.google.com/p/phyml/) set to PATH

muscle (http://www.drive5.com/muscle/) set to PATH

## Codes

### FindSeeds.py 
This code generates a seed file for TargetedAssembly.py

Useage:  python FindSeeds.py --ref_file seed.fa --ref_type nucl --fwd [options]

Options:

  -h, --help            show this help message and exit
  
  --ref_file=REF_FILE   FASTA (nucleotide or protein) file of reference sequences
  
  --ref_type=REF_TYPE   Indicate whether reference file is a nucleotide (nucl) or protein (prot) file.
  
  --fwd=FWD             FASTA/Q file of the metagenome. For paired end reads, -fwd will indicate the file with forward reads and -rev indicate the reverse reads
  
  --outFile=OUTFILE     Name of the output file that will contain the seed sequences (FASTA format).
  
  --rev=REV             FASTA/Q file that contains the reverse reads of paired end sequencing
  
  --n_threads=N_THREADS Number of threads to use [1]
  
  --evalue=EVALUE       evalue cut-off to use if --search Y [1e-5]

### TargetedAssembly.py 

A pipeline for performing targeted assembly. (Can be used without executing FindSeeds.py).

Useage: python TargetedAssembly.py --seed_file seed.fa --fwd metagenome.fa --insert_size=N --output_file=output.fasta [options]

Options:

  -h, --help                  show this help message and exit

  --seed_file=SEED_FILE       FASTA (nucleotide or protein) file of sequences to target
  
  --fwd=FWD                   FASTA/Q file of the metagenome. For paired end reads, -fwd will indicate the file with forward reads and -rev indicate the reverse reads
  
  --insert_size=INSERT_SIZE   Insert size (Required). If reads are single end (--fwd only), use the max sequence length
  
  --output_file=PRICEOUT      Name for output file. Must end in ".fasta"
  
  --rev=REV                   FASTA/Q file that contains the reverse reads of paired end sequencing
  
  --n_threads=N_THREADS       Number of threads to use [1]
  
  --mate=MATE                 Designate --mate=1 if --fwd and --rev are mate-paired (pointing away from one another) rather than paired-end (pointing towards one another)

### PredictProteins.py



### BuildTree.py


## Acknowledgements

These pipelines use software developed by others. More information is available at their resective sites (see Requirements).
