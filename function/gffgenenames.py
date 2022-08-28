# input: .gff general feature format file with gene descriptions including sequence type; enter filename in single quotes
# output: [0] array output of gene names corresponding to lines in a .gff file;
#         [1] dictionary of unique gene names and their frequencies of appearance;
#         [2] dictionary of unique tRNA anticodons and their gene copy numbers
#         [3] number of mitochondrial tRNA genes

def gffgenenames(filename):

 if "numpy" not in dir():
  import numpy

 data = numpy.genfromtxt(filename, dtype='str', delimiter='\n')
 num_lines = len(data)
 names = numpy.ndarray( num_lines, dtype='U50' )
 tRNA = numpy.array(['AAA', 'AAC', 'AAG', 'AAU', 'ACA', 'ACC', 'ACG', 'ACU', 'AGA',
                     'AGC', 'AGG', 'AGU', 'AUA', 'AUC', 'AUG', 'AUU', 'CAA', 'CAC',
                     'CAG', 'CAU', 'CCA', 'CCC', 'CCG', 'CCU', 'CGA', 'CGC', 'CGG',
                     'CGU', 'CUA', 'CUC', 'CUG', 'CUU', 'GAA', 'GAC', 'GAG', 'GAU',
                     'GCA', 'GCC', 'GCG', 'GCU', 'GGA', 'GGC', 'GGG', 'GGU', 'GUA',
                     'GUC', 'GUG', 'GUU', 'UAA', 'UAC', 'UAG', 'UAU', 'UCA', 'UCC',
                     'UCG', 'UCU', 'UGA', 'UGC', 'UGG', 'UGU', 'UUA', 'UUC', 'UUG',
                     'UUU'])
 tRNA_freqs = numpy.zeros(64, dtype='int')
 tRNA_dict = dict(zip(tRNA, tRNA_freqs))

 mito_tRNA = 0
 for i in range(0, num_lines):
  lane = data[i].split('\t')
  name = lane[2]
  names[i] = name
  if name == 'tRNA':
   if lane[0] != 'chrMito':
    anticodon_index = lane[8].find('Name=') + 8
    tRNA_type = lane[8][anticodon_index:anticodon_index + 3]
    if tRNA_type in tRNA:
     tRNA_dict[tRNA_type] = tRNA_dict[tRNA_type] + 1
   else:
    mito_tRNA = mito_tRNA + 1

 unique_names, counts_names = numpy.unique(names, return_counts=True)

 return names, dict(zip(unique_names, counts_names)), tRNA_dict, mito_tRNA
