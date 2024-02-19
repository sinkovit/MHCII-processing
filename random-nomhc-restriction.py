# Program random-nomhc-restriction.py
#
# Generate random peptides with (1) uniform amino acid usage and (2)
# amino acid usage based on frequencies in mouse proteome, and sort
# into TCR classes. Since this does not account for MHC II binding,
# this script only generates the four positions that define the TCR
# class.
#
# After are sorted into classes, construct histogram of occupancies
# (number of clases with 0, 1, 2, ... peptides/class)
#
# Author: Robert Sinkovits, San Diego Supecomputer Center

import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Generate random peptidesm without MHC II restriction')
parser.add_argument('-n', dest='npep', default=0, type=int,
                    help='Number of peptides')

args = parser.parse_args()
npep = args.npep

def fill_classes(npep, weights=None):
    tcr_class_occ = {}
    for pos1 in aa:
        for pos2 in aa:
            for pos3 in aa:
                for pos4 in aa:
                    tcr_class = pos1+pos2+pos3+pos4
                    tcr_class_occ[tcr_class] = 0

    if weights:
        fout = open("nomhc_random_weighted_classocc-hist.txt", "w")
        samp1 = np.random.choice(aa, npep, p=weights)
        samp2 = np.random.choice(aa, npep, p=weights)
        samp3 = np.random.choice(aa, npep, p=weights)
        samp4 = np.random.choice(aa, npep, p=weights)
    else:
        fout = open("nomhc_random_uniform_classocc-hist.txt", "w")
        samp1 = np.random.choice(aa, npep)
        samp2 = np.random.choice(aa, npep)
        samp3 = np.random.choice(aa, npep)
        samp4 = np.random.choice(aa, npep)

    for i in range(npep):
        samp = samp1[i] + samp2[i] + samp3[i] + samp4[i]
        tcr_class_occ[samp] += 1

    class_occ_hist = [0] * 300
    for v in tcr_class_occ.values():
        class_occ_hist[v] += 1

    for i in range(300):
        s = str(i) + " " + str(class_occ_hist[i]) + "\n"
        fout.write(s)

    fout.close()
    return
        
# Define amino acids
aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# Define mouse amino acid usage probabilities
weights = [0.068619, 0.022701, 0.048233, 0.069875, 0.036901,
           0.063779, 0.026085, 0.043845, 0.057077, 0.099819,
           0.022422, 0.035949, 0.061564, 0.047982, 0.056011,
           0.085172, 0.054170, 0.060911, 0.012034, 0.026851]

fill_classes(npep, weights=None)
fill_classes(npep, weights=weights)
