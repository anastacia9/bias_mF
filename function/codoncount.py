# input: [0] file with coding sequences alphabetically organized by gene name and in FASTA format (filename in single quotes)
# output: [0] array of arrays: the number of times each of 61 codons (stop codons not counted) appears in a gene is
#             calculated and stored in a 61-element array which in turn is stored in an n-element array (n is the number
#             of genes in the input file); array of arrays indexed by gene order of input file
#         [1] 61-element array; each element gives the probability of a codon (except stop codons) appearing in a gene
#         [2] 61-element array of codons in alphabetical order; this is the index for the 61-element arrays in [0] and in [1]

def codoncount(filename):

 if "numpy" not in dir():
  import numpy

 data = numpy.genfromtxt(filename, dtype='str', delimiter='\n')
 num_lines = len(data)
 freqs = [None] * int(num_lines/2)
 binary_vals = numpy.zeros(61)
 
 count = 0
 for j in range(1, num_lines, 2):
  cds = data[j]
  num_codons = int(len(cds)/3)
  dict_codons = {'GCT':0, 'GCC':0, 'GCA':0, 'GCG':0, 'CGT':0,
                 'CGC':0, 'CGA':0, 'CGG':0, 'AGA':0, 'AGG':0,
                 'AAT':0, 'AAC':0, 'GAT':0, 'GAC':0, 'TGT':0,
                 'TGC':0, 'CAA':0, 'CAG':0, 'GAA':0, 'GAG':0,
                 'GGT':0, 'GGC':0, 'GGA':0, 'GGG':0, 'CAT':0,
                 'CAC':0, 'ATT':0, 'ATC':0, 'ATA':0, 'TTA':0,
                 'TTG':0, 'CTT':0, 'CTC':0, 'CTA':0, 'CTG':0,
                 'AAA':0, 'AAG':0, 'ATG':0, 'TTT':0, 'TTC':0,
                 'CCT':0, 'CCC':0, 'CCA':0, 'CCG':0, 'TCT':0,
                 'TCC':0, 'TCA':0, 'TCG':0, 'AGT':0, 'AGC':0,
                 'ACT':0, 'ACC':0, 'ACA':0, 'ACG':0, 'TGG':0,
                 'TAT':0, 'TAC':0, 'GTT':0, 'GTC':0, 'GTA':0,
                 'GTG':0}
#
  for i in range(0, num_codons):
   codon = cds[3*i:3*i+3]
   if codon.find('N') == -1:   #ignore ambiguous codons
    if codon not in ['TAA', 'TAG', 'TGA']:   #stop codons
     dict_codons[codon] = dict_codons[codon] + 1   
    else:
     break
   else:
    continue
#
  key = numpy.array(list(dict_codons.keys()))
  val = numpy.array(list(dict_codons.values()))
  sort_key = key[numpy.argsort(key)]
  sort_val = val[numpy.argsort(key)]
  freqs[count] = sort_val
  binary_val = numpy.where(sort_val == 0, 0, 1)
  binary_vals = binary_vals + binary_val
  count = count + 1
#
 codon_prob = binary_vals / (num_lines/2)
 codon_freqs = numpy.array(freqs)
#
 sort_codons = numpy.array(
      ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA',
       'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC',
       'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG',
       'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT',
       'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA',
       'GTC', 'GTG', 'GTT', 'TAC', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT',
       'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT'], dtype='<U3')
#
 return codon_freqs, codon_prob, sort_codons





