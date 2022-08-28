# input: [0] file with coding sequence in FASTA format (filename in single quotes) 
#        [1] optional: codon relative adaptiveness table in dictionary format; if Sharp&Li, type 'default'
# output: [0] dictionary with gene name and corresponding CAI value;
#         [1] ndarray of gene names
#         [2] ndarray of corresponding CAI values
#         [3] number of codons factored into CAI
#         [4] dictionary with gene name and number of codons factored into CAI
#         [5] number of codons in gene
#         [6] dictionary with gene name and number of codons in gene

def CAI(filename, codontable):

 if "numpy" not in dir():
  import numpy

 data = numpy.genfromtxt(filename, dtype='str', delimiter='\n')
 
 if codontable == 'default':
  codontable = {
   'ATA':0.003, 'ATC':1.000, 'ATT':0.823, 'ATG':1.000,
   'ACA':0.012, 'ACC':1.000, 'ACG':0.006, 'ACT':0.921,
   'AAC':1.000, 'AAT':0.053, 'AAA':0.135, 'AAG':1.000,
   'AGC':0.031, 'AGT':0.021, 'AGA':1.000, 'AGG':0.003,
   'CTA':0.039, 'CTC':0.003, 'CTG':0.003, 'CTT':0.006,
   'CCA':1.000, 'CCC':0.009, 'CCG':0.002, 'CCT':0.047,
   'CAC':1.000, 'CAT':0.245, 'CAA':1.000, 'CAG':0.007,
   'CGA':0.002, 'CGC':0.002, 'CGG':0.002, 'CGT':0.137,
   'GTA':0.002, 'GTC':0.831, 'GTG':0.018, 'GTT':1.000,
   'GCA':0.015, 'GCC':0.316, 'GCG':0.001, 'GCT':1.000,
   'GAC':1.000, 'GAT':0.554, 'GAA':1.000, 'GAG':0.016,
   'GGA':0.002, 'GGC':0.020, 'GGG':0.004, 'GGT':1.000,
   'TCA':0.036, 'TCC':0.693, 'TCG':0.005, 'TCT':1.000,
   'TTC':1.000, 'TTT':0.113, 'TTA':0.117, 'TTG':1.000,
   'TAC':1.000, 'TAT':0.071, 'TAA':'---', 'TAG':'---',
   'TGC':0.077, 'TGT':1.000, 'TGA':'---', 'TGG':1.000,
  }

 num_lines = len(data)
 names = numpy.ndarray( int(num_lines/2), dtype='U15' )
 codes = numpy.zeros( int(num_lines/2) )
 length = numpy.zeros( int(num_lines/2), dtype=int )
 genelength = numpy.zeros( int(num_lines/2), dtype=int )

 for j in range(1, num_lines, 2):
  cds = data[j]
  num_codons = len(cds)/3
  w = numpy.ones(int(num_codons))
  count = 0
  n = 0
  stop = 0
  skip = 0
   
#
  for i in range(0,int(num_codons)):
   codon = cds[3*i:3*i+3]
   if codon.find('N') == -1:   #ambiguous codons
    if codon not in ['TAA', 'TAG', 'TGA']:   #stop codons
     if codon not in ['ATG', 'TGG']:   #only one codon for met and trp 
      w[i] = codontable[codon]
      count = count + 1
     else:
      skip = skip + 1
    else:
     stop = 1
     break
   else:
    n = n + 1  

#
  k = int(j/2)
  if count:
   logCAI = sum(numpy.log(w)) / count
   codes[k] = numpy.exp(logCAI)
   length[k] = count
   genelength[k] = count + n + skip + stop
  temp = str( data[j-1][1:15] )
  tloc = temp.find('\t')
  if tloc != -1:
   temp = temp[0:tloc]
  names[k] = temp
 return dict(zip(names,codes)), names, codes, length, dict(zip(names,length)), genelength, dict(zip(names,genelength))

