# Program annot_pred
#
# Description: Adds annotation to the results of MHC-II prediction output
#
# Author: Robert Sinkovits, San Diego Supecomputer Center

import argparse
import collections
import numpy as np
from collections import OrderedDict
from operator import itemgetter

parser = argparse.ArgumentParser(description='Add annotation to the MHC-II prediction output. For each accession number and entry name, append the number of peptides that bind at affinities better than the specified cutoff and the length of the protein')

parser.add_argument('-m', dest='mhcii',
                    help='Output from MHC-II binding prediction')
parser.add_argument('-c', dest='cutoff',
                    help='IC50 cutoff (nM)')
parser.add_argument('-b', dest='basename',
                    help='base name')
parser.add_argument('-f', dest='fasta', nargs='+',
                    help='fasta file')

args          = parser.parse_args()
mhcii         = args.mhcii
cutoff        = float(args.cutoff)
basename      = args.basename
fasta         = args.fasta

# Parse fasta files and create list-of-lists containing the entrynames and accession numbers
# Initialze number of binding peptides for each entryname and accession
accession_list = []
entryname_list = []
accession_cnt  = {}
entryname_cnt  = {}
accession_plen = {}
entryname_plen = {}
accession_hist = []
entryname_hist = []
index = -1

for seqfile in fasta:
    accession_list.append([])
    entryname_list.append([])
    index += 1
    with open(seqfile, 'r') as fin:
        for line in fin:
            if line.find('>') == 0:
                line = line.replace(' ', '|', 1)
                junk, accession, entryname, description = line.split('|')
                accession_list[index].append(accession)
                entryname_list[index].append(entryname)
                accession_cnt[accession]  = 0
                entryname_cnt[entryname]  = 0
                accession_plen[accession] = 0
                entryname_plen[entryname] = 0
            else:
                length = len(line.rstrip())
                accession_plen[accession] += length
                entryname_plen[entryname] += length

# Parse MHC-II prediction output and annotate with accession number
# Note that the sequence number in MHC-II prediction output must be
# decremented by one since Python arrays start with zero
index = -1
with open(mhcii, 'r') as fin:
    next(fin)
    for line in fin:
        fields = line.split()
        seq_num = int(fields[1])
        start   = int(fields[2])
        peptide  = fields[4]
        nn_core  = fields[13]
        nn_ic50  = float(fields[14])
        
        if seq_num == 1 and start == 1:
            index += 1

        if nn_ic50 <= cutoff:
            key = accession_list[index][seq_num-1]
            accession_cnt[key] += 1
            key = entryname_list[index][seq_num-1]
            entryname_cnt[key] += 1


fout = open(basename + '_accession_counts.txt', 'w')
d = OrderedDict(sorted(accession_cnt.items(), key=itemgetter(1,0)))
for key,value in d.items():
    s = key + " " + str(value) + " " + str(accession_plen[key]) + "\n"
    fout.write(s)
fout.close()

fout = open(basename + '_entryname_counts.txt', 'w')
d = OrderedDict(sorted(entryname_cnt.items(), key=itemgetter(1,0)))
for key,value in d.items():
    s = key + " " + str(value) + "\n"
    s = key + " " + str(value) + " " + str(entryname_plen[key]) + "\n"
    fout.write(s)
fout.close()


maxval = max(accession_cnt.values())
hist = [0] * (maxval + 1)
for value in accession_cnt.values():
    hist[value] += 1
fout = open(basename + '_accession_hist.txt', 'w')
for i in range(0, maxval+1):
    s = str(i) + " " + str(hist[i]) + "\n"
    fout.write(s)

maxval = max(entryname_cnt.values())
hist = [0] * (maxval + 1)
for value in entryname_cnt.values():
    hist[value] += 1
fout = open(basename + '_entryname_hist.txt', 'w')
for i in range(0, maxval+1):
    s = str(i) + " " + str(hist[i]) + "\n"
    fout.write(s)
