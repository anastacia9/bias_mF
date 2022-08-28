'n_transcripts_6179.csv' does not contain mitochondrial transcript data (genes that start with Q)
'n_transcripts_1636.csv' contains data for those genes that we also have peptide data for

The 'n' prefix signifies that the transcript data has been normalized.
 - For each strain, values were shifted up such that the minimum value is zero.
 - Transcript values were summed for each strain.
 - Each gene in each strain was divided by the sum of transcript values for that strain.
 - Each gene was multiplied by 44,155.134 (an estimate for the number of mRNA transcripts in a cell)

The code is as follows:

import pandas

df = pandas.read_csv('transcripts_6179.csv', index_col=0)
df -= df.min(axis = 0)
df = df / df.sum(axis=0)

total = df.sum(axis=0) + 44154.134

df *= total
