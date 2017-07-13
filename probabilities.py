# Program probabilities.py
#
# Append the class occupancy file (TCR_class occupancy) with the
# probabilities that a randomly selected peptide will belong to that
# class. Probabilities calculated as the product of the four amino
# acid frequencies.
#
# Author: Robert Sinkovits, San Diego Supercomputer Center

import argparse
import mhcii

parser = argparse.ArgumentParser(description='Append class probabilites, based on the product of the four amino acid frequencies, to the TCR class occupancies')

parser.add_argument('-o', dest='occupancy_file',
                    help='TCR class occupancy file')
parser.add_argument('-f', dest='aafreq_file',
                    help='Amino acid frequency file')
parser.add_argument('-r', dest='augmented_file',
                    help='TCR class occupancies augmented with probabilities')

args           = parser.parse_args()
occupancy_file = args.occupancy_file
aafreq_file    = args.aafreq_file
augmented_file = args.augmented_file

# Populate dictionary of AA frequencies
aa_freq = {}
with open(aafreq_file, 'rU') as fin:
    for line in fin:
        aa, freq = line.split()
        aa_freq[aa] = float(freq)

# Populate dictionary of TCR occupancies and calculate probabilties
tcr_occ  = {}
tcr_prob = {}
peptide_count = 0
with open(occupancy_file, 'rU') as fin:
    for line in fin:
        tcr_class, occ      = line.split()
        p1, p2, p3, p4      = list(tcr_class)
        tcr_occ[tcr_class]  = int(occ)
        tcr_prob[tcr_class] = aa_freq[p1] * aa_freq[p2] * aa_freq[p3] * aa_freq[p4]
        peptide_count += int(occ)

fout = open(augmented_file, 'w')
for k in sorted(tcr_occ.keys()):
    tcr_exp = peptide_count * tcr_prob[k]
    s = k + " " + str(tcr_occ[k]) + " " + str(tcr_exp) + "\n"
    fout.write(s)
fout.close()

print peptide_count
