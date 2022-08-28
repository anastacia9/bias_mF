# input: [0] nucleotide sequence with A's, T's, G's, and C's only;
#        [1] optional: table of codon relative adaptiveness values in dictionary format - if tAI table, type 'tAI';
# output: [0] tAI value of nucleotide sequence;

def tAI_seq(seq, codontable):

 if "numpy" not in dir():
  import numpy

 if codontable == 'tAI':
  codontable = {
  'TTT': 0.3633004926108375, 'GTT': 0.8620689655172414, 'CTT': 0.036330049261083755, 'ATT': 0.8004926108374385,
  'TGT': 0.14532019704433502, 'GGT': 0.5812807881773401, 'CGT': 0.3694581280788178, 'AGT': 0.10899014778325125,
  'TCT': 0.6773399014778326, 'GCT': 0.6773399014778326, 'CCT': 0.12315270935960593, 'ACT': 0.6773399014778326,
  'TAT': 0.29064039408867004, 'GAT': 0.5812807881773401, 'CAT': 0.2543103448275863, 'AAT': 0.3633004926108375, 
  'TTG': 0.7536945812807883, 'GTG': 0.1625615763546798, 'CTG': 0.059113300492610835, 'ATG': 0.6157635467980296,
  'TGG': 0.3694581280788178, 'GGG': 0.18226600985221678, 'CGG': 0.061576354679802964, 'AGG': 0.27832512315270935,
  'TCG': 0.1206896551724138, 'GCG': 0.09852216748768472, 'CCG': 0.19704433497536944, 'ACG': 0.14039408866995073,
  'TAG': 0.001, 'GAG': 0.39901477832512317, 'CAG': 0.23891625615763545, 'AAG': 1.0,
  'TTC': 0.6157635467980296, 'GTC': 0.6206896551724138, 'CTC': 0.061576354679802964, 'ATC': 0.5763546798029557,
  'TGC': 0.24630541871921185, 'GGC': 0.9852216748768474, 'CGC': 0.26600985221674883, 'AGC': 0.1847290640394089,
  'TCC': 0.4876847290640394, 'GCC': 0.4876847290640394, 'CCC': 0.08866995073891626, 'ACC': 0.4876847290640394,
  'TAC': 0.4926108374384237, 'GAC': 0.9852216748768474, 'CAC': 0.4310344827586207, 'AAC': 0.6157635467980296,
  'TTA': 0.4310344827586207, 'GTA': 0.12323891625615764, 'CTA': 0.1847290640394089, 'ATA': 0.12323275862068965,
  'TGA': 0.001, 'GGA': 0.1847290640394089, 'CGA': 3.6945812807877705e-05, 'AGA': 0.6773399014778326,
  'TCA': 0.18479679802955667, 'GCA': 0.3079495073891626, 'CCA': 0.6157758620689655, 'ACA': 0.24637315270935964,
  'TAA': 0.001, 'GAA': 0.8620689655172414, 'CAA': 0.5541871921182266, 'AAA': 0.4310344827586207
  }

 num_codons = int(len(seq)/3)
 w = numpy.ones(num_codons)
 count = 0
#
 n = 0 #number of ambiguous codons
 stop = 0 #number of stop codons
#
 for i in range(0, num_codons):
  codon = seq[3*i:3*i+3]
  if codon.find('N') == -1:   #ambiguous codons
   if codon not in ['TAA', 'TAG', 'TGA']:   #stop codons
    w[i] = codontable[codon]
    count = count + 1
   else:
    stop += 1
    break
  else:
   n += 1

#
# if n:
#  print(str(n) + ' ambiguous codons detected')
#
 if count:
  logtAI = sum(numpy.log(w)) / count
  mytAI = numpy.exp(logtAI)
 else:
  mytAI = None

 return mytAI
