'peptides_sum_1642.csv' contains all summed peptide data from 'peptides.txt'; peptide values associated with the same gene (in a single isolate) are summed
'peptides_sum_1637.csv' does not contain mitochondrial peptide data (genes that start with Q)
'peptides_sum_1636.csv' contains data for those genes that we also have transcript data for

The 'n' prefix signifies that the peptide data has been normalized.
 - Peptide values were summed for each strain.
 - Each gene in each strain was divided by the sum of peptide values for that strain.
 - Each gene was multiplied by 53e6 (an estimate for the number of proteins in a cell)

See README.txt in yeast_transcripts folder for code.
