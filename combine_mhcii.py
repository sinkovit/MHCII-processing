# Program combine_mhcii.py
#
# Description: Combine the results from a pair of MHC Class II peptide
# binding prediction analyses.
#
# Input:
# (1) Results of MHC class II prediction for foreign proteome
# (2) Results of MHC class II prediction for self proteome
# (3) IC50 cutoff
# (4) Method (NN or SMM)
#
# Output:
#
# (1) MHC class II binding prediction results for foreign proteome
# annotated by TCR class and occupancy from reference proteome. The
# method and IC50 cutoff are used both to filter the results from the
# foreign proteome and to determine the class occupancy for the
# reference proteome. Note that both self and foreign predictions need
# to be based on same MHC II molecule

import argparse
import mhcii

parser = argparse.ArgumentParser(description='Annotate results of foreign MHC II peptide binding prediction with TCR occupancies from self binding predictions. Done using a specified method and IC50 cutoff')

parser.add_argument('-s', dest='self_prediction',
                    help='MHC II predictions for self peptides')
parser.add_argument('-f', dest='foreign_prediction',
                    help='MHC II predictions for foreign peptides')
parser.add_argument('-r', dest='results',
                    help='Results for foreign peptide binding relative to self')
parser.add_argument('-m', dest='method', default='nn',
                    help='MHC II prediction method (nn or smm)')
parser.add_argument('-c', dest='cutoff', type=int,
                    help='Cutoff affinity (IC50) in nM')

args                = parser.parse_args()
self_prediction     = args.self_prediction
foreign_prediction  = args.foreign_prediction
results             = args.results
method              = args.method
cutoff              = args.cutoff

print 'Self proteome prediction file:    ', self_prediction
print 'Foreign proteome prediction file: ', foreign_prediction
print 'MHC binding prediction method:    ', method
print 'IC50 cuotff:                      ', cutoff

# Determine TCR class occupancies for self peptides
print
print '-- Self-binding statistics --'
occupancy = {}
if method == 'smm':
    nn_file, smm_file = mhcii.mhcii_to_TCR_classes(self_prediction, 0, cutoff)
    self_occ = smm_file
else:
    nn_file, smm_file = mhcii.mhcii_to_TCR_classes(self_prediction, cutoff, 0)
    self_occ = nn_file

with open(self_occ, 'r') as fin:
    for line in fin:
        tcr_class, tcr_occ = line.split()
        occupancy[tcr_class] = tcr_occ


# Process the foreign peptide binding file and write annotated results
counter_foreign = 0
counter_foreign_binding = 0

fout = open(results, 'w')
with open(foreign_prediction, 'r') as fin:
    next(fin)
    for line in fin:
        counter_foreign += 1
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

        if method == 'smm' and smm_ic50 < cutoff:
            counter_foreign_binding += 1
            s = peptide + " " + smm_core + " " + tcr_smm_class + " " + str(smm_ic50) + " " + occupancy[tcr_smm_class] + "\n"
            fout.write(s)
        if method == 'nn' and nn_ic50 < cutoff:
            counter_foreign_binding += 1
            s = peptide + " " +  nn_core + " " + tcr_nn_class + " "  + str(nn_ic50)  + " " + occupancy[tcr_nn_class]  + "\n"
            fout.write(s)

fout.close()

pct_binding = 100 * float(counter_foreign_binding)/float(counter_foreign)
print
print '-- Foreign-binding statistics --'
print
print '# peptides:               ', counter_foreign
print '# peptides IC50 < cutoff: ', counter_foreign_binding
print 'Percent binding:          ', pct_binding

