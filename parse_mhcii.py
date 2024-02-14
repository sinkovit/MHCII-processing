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


def mhcii_to_TCR_classes(infile, nn_cutoff, smm_cutoff):
    # Initialize TCR classes
    smm_dict = {}
    nn_dict  = {}
    aas = 'ACDEFGHIKLMNPQRSTVWY'
    for pos1 in aas:
        for pos2 in aas:
            for pos3 in aas:
                for pos4 in aas:
                    tcr_class = pos1+pos2+pos3+pos4
                    smm_dict[tcr_class] = 0
                    nn_dict[tcr_class]  = 0


    # Initialize counts for peptides with bindings below affinity
    # threshold, total number of peptides and number of empty classes

    nn_count_below_IC50  = 0
    smm_count_below_IC50 = 0
    smm_empty_classes    = 0
    nn_empty_classes     = 0
    peptide_count        = 0

    # Open and parse output from MHC II prediction program
    with open(infile, 'r') as fin:
        for line in fin:
            # Skipping header lines
            if line.find('allele') >= 0:
                continue

            peptide_count += 1
            fields = line.split()

            # Get information on peptide, SMM scoring and NN scoring
            peptide  = fields[4]
            smm_core = fields[10]
            smm_ic50 = float(fields[11])
            smm_rank = float(fields[12])
            nn_core  = fields[13]
            nn_ic50  = float(fields[14])
            nn_rank  = float(fields[15])

            # Determine the TCR class from core nonamer positions 2,3,5 and 8
            # Note that code reflects Python counting from zero
            tcr_smm_class = smm_core[1]+smm_core[2]+smm_core[4]+smm_core[7]
            tcr_nn_class  = nn_core[1]+nn_core[2]+nn_core[4]+nn_core[7]

            if smm_ic50 < smm_cutoff:
                smm_dict[tcr_smm_class] += 1
                smm_count_below_IC50 += 1
            if nn_ic50 < nn_cutoff:
                nn_count_below_IC50 += 1
                nn_dict[tcr_nn_class]   += 1

    # Write results for SMM method
    smm_file = 'smm_'+str(smm_cutoff)+'.out'
    if smm_cutoff > 0:
        sorted_keys = sorted(smm_dict.keys())
        smm_empty_classes = 0
        fout = open(smm_file, 'w')
        for k in sorted_keys:
            if smm_dict[k] == 0:
                smm_empty_classes += 1
            s = k + " " + str(smm_dict[k]) + "\n"
            fout.write(s)
        fout.close()

        pct_binding = 100 * float(smm_count_below_IC50)/float(peptide_count)
        print()
        print('SMM statistics')
        print('IC50 cutoff:              ', smm_cutoff)
        print('# peptides:               ', peptide_count)
        print('# peptides IC50 < cutoff: ', smm_count_below_IC50)
        print('Percent binding:          ', pct_binding)

# Write results for NN method
    nn_file = 'nn_'+str(nn_cutoff)+'.out'
    if nn_cutoff > 0:
        sorted_keys = sorted(nn_dict.keys())
        nn_empty_classes = 0
        fout = open(nn_file, 'w')
        for k in sorted_keys:
            if nn_dict[k] == 0:
                nn_empty_classes += 1
            s = k + " " + str(nn_dict[k]) + "\n"
            fout.write(s)
        fout.close()

        pct_binding = 100 * float(nn_count_below_IC50)/float(peptide_count)
        print()
        print('NN statistics')
        print('IC50 cutoff:              ', nn_cutoff)
        print('# peptides:               ', peptide_count)
        print('# peptides IC50 < cutoff: ', nn_count_below_IC50)
        print('Percent binding:          ', pct_binding)

    return nn_file, smm_file

mhcii_to_TCR_classes(infile, nn_cutoff, smm_cutoff)
