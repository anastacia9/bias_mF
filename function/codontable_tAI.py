# input: .gff general feature format file with gene descriptions including sequence type; enter filename in single quotes
# output: codon table in dictionary format for tAI calculations
#
# Note: gffgenenames.py needs to be in working directory

def codontable_tAI(filename):

 if "gffgenenames" not in dir():
  import gffgenenames

 if "numpy" not in dir():
  import numpy

 names = gffgenenames.gffgenenames(filename)
 tRNA = numpy.array(list(names[2].keys()))
 rev_translate = {'A':'T', 'U':'A', 'C':'G', 'G':'C'}

 codon = numpy.array([rev_translate[tRNA[i][2]] + rev_translate[tRNA[i][1]] + rev_translate[tRNA[i][0]] for i in range(0, 64)])
 W = numpy.zeros(64)
 
 s_UA = 0
 s_UG = 0.41
 s_CG = 0
 s_CA = 0.28
 s_AU = 0
 s_AA = 0.9999
 s_GC = 0
 s_GU = 0.68

 for k in range(0, 64):
  first = tRNA[k][0]
  second = tRNA[k][1]
  third = tRNA[k][2]
  if first == 'A':
   W[k] = (1-s_UA)*names[2][tRNA[k]] + (1-s_UG)*names[2]['G'+second+third]
  if first == 'G':
   W[k] = (1-s_CG)*names[2][tRNA[k]] + (1-s_CA)*names[2]['A'+second+third]
  if first == 'U':
   W[k] = (1-s_AU)*names[2][tRNA[k]] + (1-s_AA)*names[2]['A'+second+third]
  if first == 'C':
   W[k] = (1-s_GC)*names[2][tRNA[k]] + (1-s_GU)*names[2]['U'+second+third]

 index = numpy.where(codon == 'ATG')[0][0]
 W[index] = (1-s_GC)*names[2][tRNA[index]]

 w = W / max(W)

 return dict(zip(codon, w)), dict(zip(codon, W))
