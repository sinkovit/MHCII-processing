# Program aa_stats.py
#
# Description: Parses FASTA file containing protein sequences to
# determine amino acid counts and probabilities. Typically used for
# analysis of full proteomes. Produces two sets of outputs based on
#
# (1) 20 standard residues
# (2) 20 standard + 6 ambiguous residues
#
# Author: Robert Sinkovits, San Diego Supecomputer Center

import argparse

parser = argparse.ArgumentParser(description='Determine amino acid frequencies in FASTA file')

parser.add_argument(dest='infile',
                    help='Input file in FASTA format')

args          = parser.parse_args()
infile        = args.infile

# Define and initialize dictionary of AAs
aa_count = {}
std_aas = 'ACDEFGHIKLMNPQRSTVWY'
all_aas = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
for c in all_aas:
    aa_count[c] = 0

# Loop over records in proteome and accumulate AA counts
fin = open(infile, 'r')
for line in fin:
    if line[0] == '>':
        continue
    for c in list(line.strip('\n')):
        aa_count[c] += 1
fin.close()

# Get count for standard AAs and all AAs
count_std = 0
for c in std_aas:
    count_std += aa_count[c]

count_all = 0
for c in all_aas:
    count_all += aa_count[c]

print('# residues:     ', count_all)
print('# std residues: ', count_std)
print()

# Print out results for standard AAs
print("Statistics for standard AAs")
print("residue  count       frequency")
for c in std_aas:
    freq = float(aa_count[c]) / float(count_std)
    print( '{}     {:10d}    {:9.6f}'.format(c, aa_count[c], freq))

# Print out results for standard + ambiguous AAs
print()
print("Statistics for standard + ambiguous AAs")
print("residue  count       frequency")
for c in all_aas:
    freq = float(aa_count[c]) / float(count_all)
    print('{}     {:10d}    {:9.6f}'.format(c, aa_count[c], freq))

