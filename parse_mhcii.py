# Program parse_mhcii.py
#
# Description: Parse output from MHC Class II peptide binding
# prediction program (mhc_II_binding.py) and sort peptides with
# affinities below a specified threshold into T cell receptor
# classes. Note that this does not work with the consensus method
# since different algorithms may predict different core nonamers for a
# given 15-mer.
#
# Author: Robert Sinkovits, San Diego Supercomputer Center

import argparse
import mhcii

parser = argparse.ArgumentParser(description='Parse MHC II prediction data and sort peptides into T cell receptor classes')
parser.add_argument(dest='infile',
                    help='Input file in FASTA format')
parser.add_argument('-n', dest='nn_cutoff', default=0, type=int,
                    help='Cutoff affinity for NN algorithm affinity prediction')
parser.add_argument('-s', dest='smm_cutoff', default=0, type=int,
                    help='Cutoff affinity for SMM algorithm affinity prediction')

args          = parser.parse_args()
infile        = args.infile
nn_cutoff     = args.nn_cutoff
smm_cutoff    = args.smm_cutoff

mhcii.mhcii_to_TCR_classes(infile, nn_cutoff, smm_cutoff)
