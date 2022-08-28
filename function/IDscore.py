# input: [0] unwrapped MUSCLE alignments file in .fa or .afa format; assumed equal number of sequences aligned for each gene; limit is 1e4 genes as of 08/03/18;
#        [1] number of sequences aligned for each gene
# output: [0] vector of gene names;
#         [1] vector of percent identity scores

def IDscore(filename, num):

 if "numpy" not in dir():
  import numpy

 data = numpy.genfromtxt(filename, dtype='str', delimiter='\n')

 scores = numpy.zeros(int(len(data)/(2*num)))
 names = numpy.ndarray(int(len(data)/(2*num)), dtype='U20')

# Compute alignment scores for each gene
 count = 0
 for i in range (1, len(data), 2*num):
  # Set up matrix for all sequences corresponding to gene
  length = len(data[i])
  mat = numpy.ndarray(shape = [num, length], dtype = 'int')
  for j in range(0, num):
   # Fill matrix rows with characters of each sequence
   sequence = numpy.array( [data[i+2*j][k] for k in range(0, length)], dtype='<U1' )
   numeric = numpy.fromstring(sequence, dtype = 'int')
   mat[j] = numeric
  # Populate scores and names vector with data
  temp = numpy.equal( numpy.amin(mat, axis=0), numpy.amax(mat, axis=0) )
  scores[count] = sum(temp) / len(temp) * 100
  name = data[i-1][1:21]
  if "trans_" in name:
   name = name[6:]
  index = name.find(' ')
  if index != -1:
    name = name[:index]
  names[count] = name
  count = count + 1
 
 return names, scores


### REFERENCE: numbers corresponding to letters of alphabet
# a = numpy.array(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'])
# numeric = numpy.fromstring(a, dtype = 'int')
# print(numeric)
# print(numeric)
# [ 97  98  99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122]
