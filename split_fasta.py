# Program split_fasta.py
#
# Description: Standalone Python program to split a FASTA file into
# chunks of a specified size. Also has capabilites to 
#
# (1) Omit sequences below a given threshold length
# (2) Remove '>' embedded within the annotation line
# (3) Trim ambiguous residues from front and/or end of sequence
# (4) Omit sequences that contain ambiguous residues (after trimming)
#
# Author: Robert Sinkovits, San Diego Supercomputer Center

import re
import argparse

def write_fasta_record(sequence, header, minlen, trim_nonstd, 
                       reject_nonstd, nonstd_chars, fout, ferr, 
                       nshort, nnonstd):
    if (len(sequence) > 0):
        # Remove trailing whitespace
        sequence = sequence.rstrip()
        # Remove leading and trailing nonstandard characters
        if (trim_nonstd):
            if any(c in sequence for c in nonstd_chars):
                sequence = sequence.rstrip(nonstd_chars)
            if any(c in sequence for c in nonstd_chars):
                sequence = sequence.lstrip(nonstd_chars)
        # Replace newline
        sequence += '\n'

        # If sequence is too short or contains non-standard
        # residues, write to error file. Otherwise write to
        # output file
        if len(sequence) < minlen + 1:
            ferr.write(header)
            ferr.write(sequence)
            nshort += 1
        elif any(c in sequence for c in nonstd_chars) and reject_nonstd:
            ferr.write(header)
            ferr.write(sequence)
            nnonstd += 1
        else:
            fout.write(header)
            fout.write(sequence)

    return nshort, nnonstd

# Command line argument parsing
parser = argparse.ArgumentParser(description='Split FASTA file with optional additional processing')

parser.add_argument(dest='infile',
                    help='Input file in FASTA format')
parser.add_argument('-c', dest='chunk', default=10000000000, type=int,
                    help='Specify number of sequences assigned to each file')
parser.add_argument('-m', dest='minlen', default=0, type=int,
                    help='Omit sequences shorter than MINLEN (after trimming)')
parser.add_argument('-r', dest='rmgt', action='store_true',
                    help='Remove embedded ">" from annotation records')
parser.add_argument('-t', dest='trim_nonstd', action='store_true',
                    help='Trim non-std characters from start/end of sequence')
parser.add_argument('-x', dest='reject_nonstd', action='store_true',
                    help='Reject sequences w/ non-std characters (after trimming)')
parser.add_argument('-N', dest='nucleotide', action='store_true',
                    help='Input file contains nucleotide sequences')
parser.add_argument('-A', dest='amino_acid', action='store_true',
                    help='Input file contains mino acid sequences (default)')

args          = parser.parse_args()
infile        = args.infile
chunk         = args.chunk
rmgt          = args.rmgt
minlen        = args.minlen
trim_nonstd   = args.trim_nonstd
reject_nonstd = args.reject_nonstd
nucleotide    = args.nucleotide
amino_acid    = args.amino_acid

# Define nonstandard residues and bases
if nucleotide:
    nonstd_chars = 'bdhkmnrsvwxyBDHKMNRSVWXY-'
else:
    nonstd_chars =  'BJOUXZ'

# Initializations
seqfile_suffix = 0
seqfile_name = 'seq_' + '%03d'%seqfile_suffix + '.fasta'
fout = open(seqfile_name, 'w')
ferr = open('errors.fasta', 'w')
count = -1
header = ''
sequence = ''
nshort = 0
nnonstd = 0

# Loop over records in FASTA file
fin = open(infile, 'rU')
for line in fin:
    if line[0] == '>':
        count += 1

        nshort, nnonstd = write_fasta_record(sequence, header, minlen, trim_nonstd, 
                                             reject_nonstd, nonstd_chars, fout, ferr, 
                                             nshort, nnonstd)

        # Remove embedded '>' in annotation line
        if (rmgt):
            line = re.sub('>', '', line)
            line = '>' + line

        # Store the annotation
        header = line

        # Reinitialize the sequence to empty string
        sequence = ''
    else:
        sequence = sequence + line
    
    # Close current file and open new file for next chunk
    if (count == chunk):
        fout.close()
        seqfile_suffix += 1
        seqfile_name = 'seq_' + '%03d'%seqfile_suffix + '.fasta'
        fout = open(seqfile_name, 'w')
        count = 0
        

# Write out the last FASTA record
nshort, nnonstd = write_fasta_record(sequence, header, minlen, trim_nonstd, 
                                     reject_nonstd, nonstd_chars, fout, ferr, 
                                     nshort, nnonstd)
        

# Write out some statistics
print '# reject sequences below minimum length:   ', nshort
print '# rejected sequences w/ nonstd characters: ', nnonstd
print
