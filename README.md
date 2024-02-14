# MHC-II-processing

This repository contains a collection of scripts and utilities that
were used for generating data in the publication *PD1 and CD73
synergistically limit CD4 T cell responses to autoantigens*

### Step 1 - Preprocess and optionally split proteome into chunks

The MHC II binding prediction tool requires that the input sequences
contain at least 15 residues. In addition, the sequences cannot
contain ambiguous residues (BJOXZ) or selenocysteine (U), nor can the
annotation lines contain embedded '>' characters.

split-and-preprocess.py addresses all of these issues in addition to
optionally splitting the input file into chunks. The following command
was used for preprocessing the mouse proteome and generating chunks of

2000 records; sequences that do not satisfy length and residue
requirements are written to errors.fasta. Run with -h to see usage.

```
python3 split-and-preprocess.py -c 2000 -m 15 -r -t -x -A uniprot-proteome%3AUP000000589.fasta
```

### Step 2 - Run MHC II binding prediction tool

Run binding prediction tool on each chunk of the mouse proteome file
using MHC II allele $I-a^b$. The example below does this for a small
sample (first three records) of the first chunk generated in the
previous step. Syntax is for version 2.13 of IEDB tools. 

```
python2 [path to mhc_ii installation]/mhc_II_binding.py IEDB_recommended H2-IAb seq_000.fasta > seq_000.txt
```
