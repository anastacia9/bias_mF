# Function for calculating the percent GC content of an input cDNA sequence.
#
# input:  [0] sequence (string of A's, T's, C's, and G's)
# output: [0] % GC content

def GC(sequence):
 G = sequence.count('G')
 C = sequence.count('C')
 A = sequence.count('A')
 T = sequence.count('T')

 if G+C+A+T != len(sequence):
  print('non G, C, A, T bases in sequence')

 percent = (G+C) / len(sequence) * 100

 return(percent)
