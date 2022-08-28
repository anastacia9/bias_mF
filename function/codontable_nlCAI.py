# Compute a median table of nlCAI codon relative adaptiveness values from genomic sequence files in input path.
# This algorithm calculates a codon i's relative adaptiveness value by (1) calculating the fraction of codons that are
# codon i in each reference set gene, (2) adding up all such fractions across reference set genes, and (3) dividing this
# sum of fractions for codon i by the sum of fractions for the most prevalent synonymous codon.
#
# input:  [0] path to unwrapped genomic data files in FASTA format; order of genes must be identical;
#             destination must contain necessary genomic data files only;
#         [1] numpy array of reference set genes
#
# output: [0] dictionary of median codon relative adaptiveness values across strains
#
# Note: Carbone array of reference set genes: carbone_refset = numpy.array(['YJL052W', 'YJR009C', 'YGR192C', 'YLR044C', 'YGR254W', 'YHR174W', 'YKL060C', 'YDR050C', 'YKL152C', 'YOL086C', 'YLR355C', 'YAL038W', 'YBR118W', 'YPR080W', 'YOR133W', 'YDR385W', 'YLR249W', 'YEL034W', 'YDL229W', 'YLL024C', 'YCR012W', 'YLR110C', 'YMR116C', 'YJL189W', 'YGL030W', 'YLR061W', 'YBR181C', 'YOR369C', 'YPR043W', 'YLL045C', 'YPL131W', 'YJR123W', 'YHL033C', 'YGL135W', 'YPL090C', 'YDR012W', 'YLR075W', 'YOR063W', 'YPL220W', 'YHL015W', 'YOL120C', 'YER074W', 'YKL180W', 'YBR031W', 'YJL190C', 'YOR293W', 'YCR031C', 'YBR189W', 'YBL092W', 'YLR340W', 'YOL039W', 'YNL178W', 'YML024W', 'YGL123W', 'YPL249C-A', 'YOL121C', 'YPR132W', 'YLR167W', 'YLR029C', 'YGL189C', 'YFR031C-A'])

def codontable_nlCAI(path, array):

 if "os" not in dir():
  import os

 if "numpy" not in dir():
  import numpy

 # Define amino acids and their synonymous codons.
 Ala = numpy.array([ 'GCT', 'GCC', 'GCA', 'GCG' ])
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
 amino_flatten = numpy.hstack(amino)
 dict_rav = {}

 # Find index of each reference set gene name in files.
 num_refset_genes = len(array)
 file = os.listdir(path)[0]
 data = numpy.genfromtxt(path + '\\' + file, dtype='str', delimiter='\n')
 names_indices = numpy.zeros(num_refset_genes, dtype='int')
 for i in range(0, num_refset_genes):
  names_indices[i] = numpy.where(data == '>' + array[i])[0][0]

 # Find index of each reference set gene sequence in first file.
 seqs_indices = names_indices + 1

 # Specify empty array of arrays of strain codon fraction sums for median-ing later on.
 num_files = len(os.listdir(path))
 dict_w_aggregate = [None] * num_files

 # One file at a time . . .
 for file in range(0, num_files):
  # compute reference set codon fractions . . .
  codon_fracs = numpy.zeros(61)
  data = numpy.genfromtxt(path + '\\' + os.listdir(path)[file], dtype='str', delimiter='\n')
  seqs = data[seqs_indices]
  # for each gene sequence . . .
  for seq in seqs:
   num_codons = int(len(seq)/3)
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
   # by populating a codon frequencies dictionary, . . . 
   for k in range(0, num_codons):
    codon = seq[3*k:3*k+3]
    if codon.find('N') == -1:
     if codon not in ['TAA', 'TAG', 'TGA']:
      dict_codons[codon] = dict_codons[codon] + 1
     else:
      break
    else:
     continue
 
   # dividing the dictionary by the total number of codons in the sequence, . . .
   seq_codon_fracs = numpy.array(list(dict_codons.values()))/sum( numpy.array(list(dict_codons.values())) )
 
   # adding sequence codon fractions together.
   codon_fracs = codon_fracs + seq_codon_fracs

  # Divide sum of codon i fractions by sum of fractions of most prevalent synonymous codon.
  dict_codon_fracs = dict(zip(amino_flatten, codon_fracs))
  dict_w = {}
  for i in amino:
   frequency = numpy.asarray([dict_codon_fracs[x] for x in i])
   w = frequency / max(frequency)
   dict_w = dict(dict_w, **dict(zip(i, w)))
   
  dict_w_aggregate[file] = numpy.array(list( dict_w.values() ))

 dict_rav = dict(zip(  numpy.array(list(dict_w.keys())), numpy.median(dict_w_aggregate, axis=0)  ))

 for key, value in dict_rav.items():
  if value < 0.001:
   dict_rav[key] = 0.001

 return dict_rav
