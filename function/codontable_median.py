# Function for calculating a table of codon values for each genomic sequence file in input path.
# Computations are based on a given amount of highly expressed genes.
# input: [0] path to unwrapped genomic data files saved in FASTA format; genes in same order as in csv file
#        [1] number of most highly expressed genes to consider as reference set
#        [2] csv file with expression data; first column has gene name labels and successive columns are labelled by strain 
#            or file name and contain expression data for the corresponding file in directory; (csv file in single quotes)
#            
# output: [0] dictionary of median relative adaptiveness values (each codon has a median value across strains)
#         [1] dictionary of relative adaptiveness value tables (dictionary keys are strain names; dictionary values are 
#             the corresponding relative adaptiveness tables)
#         [2] list of gene names used for calculating table values
#

def codontable_median(path, num_genes, csvfile):


# Define numpy array of synonymous codons for each of 20 amino acids.

 if "numpy" not in dir():
  import numpy

 Ala = numpy.array([ 'GCT', 'GCC', 'GCA', 'GCG'  ])
 Arg = numpy.array([ 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG' ])
 Asn = numpy.array([ 'AAT', 'AAC' ])
 Asp = numpy.array([ 'GAT', 'GAC' ])
 Cys = numpy.array([ 'TGT', 'TGC' ])
 Gln = numpy.array([ 'CAA', 'CAG' ])
 Glu = numpy.array([ 'GAA', 'GAG' ])
 Gly = numpy.array([ 'GGT', 'GGC', 'GGA', 'GGG' ])
 His = numpy.array([ 'CAT', 'CAC' ])
 Ile = numpy.array([ 'ATT', 'ATC', 'ATA' ])
 Leu = numpy.array([ 'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG' ])
 Lys = numpy.array([ 'AAA', 'AAG' ])
 Met = numpy.array([ 'ATG' ])
 Phe = numpy.array([ 'TTT', 'TTC' ])
 Pro = numpy.array([ 'CCT', 'CCC', 'CCA', 'CCG' ])
 Ser = numpy.array([ 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC' ])
 Thr = numpy.array([ 'ACT', 'ACC', 'ACA', 'ACG' ])
 Trp = numpy.array([ 'TGG' ])
 Tyr = numpy.array([ 'TAT', 'TAC' ])
 Val = numpy.array([ 'GTT', 'GTC', 'GTA', 'GTG' ])
 amino = numpy.array([Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile, Leu, Lys, Met, Phe, Pro, Ser, Thr, Trp, Tyr, Val])
 dict_rav = {}

# Rank genes by median expression level.

 if "pandas" not in dir():
  import pandas

 expression = pandas.read_csv(csvfile, index_col=0)
 list_medis = numpy.array(list(expression.median(axis=1)))
 list_genes = numpy.array(list(expression.index))
 index = numpy.argsort(list_medis)
 sort_list_medis = numpy.array(list_medis[index])
 sort_list_genes = numpy.array(list_genes[index])
 medis = sort_list_medis[-num_genes:]
 genes = sort_list_genes[-num_genes:]

# Load data and start for-loop.

 if "os" not in dir():
  import os

 directory = os.fsencode(path)
 count = 0
 for file in os.listdir(directory):
  filename = os.fsdecode(file)
  data = numpy.genfromtxt(path + '\\' + filename, dtype='str', delimiter='\n')
  names = data[0::2]
  index0 = [2*numpy.where(names == '>' + i)[0][0]+1 for i in genes]
  seqs = data[index0]

# Count up how many times each of 61 codons appears in the input file (disregard stop codons). Store result in dictionary format.
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
  for j in range(0, num_genes):
   cds = seqs[j]
   for i in range(0, int( len(cds)/3 )):
    codon = cds[3*i:3*i+3]
    if codon.find('N') == -1:
     if codon.find('TAA') == -1:
      if codon.find('TAG') == -1:
       if codon.find('TGA') == -1:
        if codon in dict_codons:
         dict_codons[codon] = dict_codons[codon] + 1
        else:
         print('unexpected codon at gene index: ' + j + ' and codon index: ' + i + ' in file: ' + filename)
       else:
        break
      else:
       break 
     else:
      break


# Calculate relative adaptiveness values.
# Relative adaptiveness value = (number of appearances of synonymous codon) / (number of appearances of most common synonymous codon)
  dict_w = {}
  for i in amino:
   frequency = numpy.asarray([ dict_codons[x] for x in i ])
   w = frequency / max(frequency)
   dict_w = dict(dict_w, **dict(zip(i, w)))
   dict_rav[filename[:-11]] = dict_w 
  count = count + 1
 
# Calculate the median table of codon values.
 dict_med = {}
 for codon in numpy.concatenate(amino).ravel():
  dict_med[codon] = numpy.median([dict_rav[strain][codon] for strain in expression.columns.values])
 
 return dict_med, dict_rav, genes
