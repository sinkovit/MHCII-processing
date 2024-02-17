# ------------------------------------------------------------------------
# The following steps are computationally intensive and we recommend
# using data downloaded from FigShare DOI 10.6084/m9.figshare.15057975
# ------------------------------------------------------------------------

# Preprocess mouse proteome, split into chunks, and run the MHC II
# binding prediction program on each chunk. After calculations finish,
# combine results into a single file. With '-c 2000' this will result
# in 26 files of 2000 sequence each. Adjust value to increase or
# decrease number of files.

###python3 split-and-preprocess.py -c 2000 -m 15 -r -t -x -A uniprot-proteome.3AUP000000589.fasta
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_000.fasta > seq_000.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_001.fasta > seq_001.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_002.fasta > seq_002.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_003.fasta > seq_003.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_004.fasta > seq_004.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_005.fasta > seq_005.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_006.fasta > seq_006.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_007.fasta > seq_007.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_008.fasta > seq_008.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_009.fasta > seq_009.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_010.fasta > seq_010.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_011.fasta > seq_011.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_012.fasta > seq_012.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_013.fasta > seq_013.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_014.fasta > seq_014.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_015.fasta > seq_015.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_016.fasta > seq_016.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_017.fasta > seq_017.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_018.fasta > seq_018.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_019.fasta > seq_019.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_020.fasta > seq_020.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_021.fasta > seq_021.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_022.fasta > seq_022.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_023.fasta > seq_023.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_024.fasta > seq_024.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_025.fasta > seq_025.txt
###cat seq_???.txt > mhciipred.mouse.UP000000589.txt

# Repeat for synthetic mouse genome with same number of sequences and
# sequence lengths, except where residues are randomly selected based
# on the mouse amino acide frequencies.

###python3 split-and-preprocess.py -c 2000 -m 15 -r -t -x -A uniprot-proteome_random.3AUP000000589.fasta
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_000.fasta > seq_000.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_001.fasta > seq_001.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_002.fasta > seq_002.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_003.fasta > seq_003.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_004.fasta > seq_004.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_005.fasta > seq_005.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_006.fasta > seq_006.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_007.fasta > seq_007.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_008.fasta > seq_008.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_009.fasta > seq_009.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_010.fasta > seq_010.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_011.fasta > seq_011.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_012.fasta > seq_012.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_013.fasta > seq_013.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_014.fasta > seq_014.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_015.fasta > seq_015.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_016.fasta > seq_016.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_017.fasta > seq_017.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_018.fasta > seq_018.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_019.fasta > seq_019.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_020.fasta > seq_020.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_021.fasta > seq_021.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_022.fasta > seq_022.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_023.fasta > seq_023.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_024.fasta > seq_024.txt
###python2 mhc_II_binding.py IEDB_recommended H2-IAb seq_025.fasta > seq_025.txt
###cat seq_???.txt > mhciipred.mouse_random.UP000000589.txt



# ------------------------------------------------------------------------
# The previous steps are computationally intensive and we recommend
# using data downloaded from FigShare DOI 10.6084/m9.figshare.15057975
# ------------------------------------------------------------------------

# Calculate class occupancies and histogram for mouse proteome
python3 parse_mhcii.py -n 1000 -p mouse mhciipred.mouse.UP000000589.txt

# Calculate class occupancies and histogram for synthetic mouse proteome
python3 parse_mhcii.py -n 1000 -p mouse_random mhciipred.mouse_random.UP000000589.txt

# Calculate amino acid usage for mouse proteome
python3 aa_stats.py uniprot-proteome%3AUP000000589.fasta
