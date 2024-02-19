# MHC-II-processing

This repository contains a collection of scripts and utilities that
were used for generating data in the publication *PD1 and CD73
synergistically limit CD4 T cell responses to autoantigens*

The full output from the MHC II binding predictions, especially for
mouse self-peptide binding, is too large to store in GitHub. The
examples below involve very small subsets of the data. The complete
data, including the version of the mouse proteome used in this study
are published on FigShare at the DOI 10.6084/m9.figshare.15057975

The workflow using the complete data is given in MHCII-workflow.sh

### Step 1 - Preprocess and optionally split proteome into chunks

The MHC II binding prediction tool requires that the input sequences
contain at least 15 residues. In addition, the sequences cannot
contain ambiguous residues (BJOXZ) or selenocysteine (U), nor can the
annotation lines contain embedded '>' characters.

`split-and-preprocess.py` addresses all of these issues in addition to
optionally splitting the input file into chunks. The following command
was used for preprocessing the mouse proteome and generating chunks of
2000 records; sequences that do not satisfy length and residue
requirements are written to `errors.fasta`. Run with `-h` to see usage.

```
python3 split-and-preprocess.py -c 2000 -m 15 -r -t -x -A uniprot-proteome%3AUP000000589.fasta
```

To demonstrate how `split-and-preprocess.py` handles short lines,
ambiguous characters, etc., execute the following

```
python3 split-and-preprocess.py -m 15 -r -t -x -A data/aa.fasta
python3 split-and-preprocess.py -m 15 -r -t -x -N data/bases.fasta
```

### Step 2 - Run MHC II binding prediction tool

Run binding prediction tool on each chunk of the mouse proteome file
using MHC II allele $IA^b$. The example below does this for a small
sample (first three records) of the first chunk generated in the
previous step. Syntax is for version 2.13 of IEDB tools. 

```
python2 [path to mhc_ii installation]/mhc_II_binding.py IEDB_recommended H2-IAb data/seq_000.fasta > data/seq_000.txt
```

### Step 3 - Calculate self-peptide class occupancies and histograms

After running the MHC II binding prediction tool, concatenate the
output and use `parse_mhcii.py` to calculate the self-peptide class
occupancies based on the NN algorithm. The script automatically strips
out multiple header lines and the affinity threshold is in units of
nM. Run with `-h` to see usage. In the example below, output files
will be named `mouse_nn_1000nM_classocc.txt` and
`mouse_nn_1000nM_classocc-hist.txt`.

Note that in this small example, most classes will have zero
occupancy. To see occupied classes, run `grep -v 0 mouse_nn_1000nM_classocc.txt`

```
python3 parse_mhcii.py -n 1000 -p mouse data/seq_000.txt
```

### Step 4 - Calculate amino acid usage stats for mouse proteome

`aa_stats.py` gives amino acid usages for both standard and
non-standard (BOUJXZ) residues across entire proteome

```
python3 aa_stats.py uniprot-proteome%3AUP000000589.fasta
```

### Step 5 - Generate histograms for random peptides w/o MHC restriction

```python3 random-nomhc-restriction.py -n 2516523```
