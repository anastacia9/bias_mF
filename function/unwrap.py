# input: file in fasta format; filename in single quotes
# output: 'unwrap.afa' file with gene name and corresponding sequence on alternating lines
# note: longest gene in the yeast genome is MDN1 (YLR106C) which is ~14.7 kb long; if sequence of gene longer than 15 kb, unwrap.py fails
# note: if name of gene longer than 40 characters, gene names in unwrap.afa will be different than in input file

def unwrap(filename):

 if "numpy" not in dir():
  import numpy

 names = numpy.ndarray( int(4e4), dtype='U60' )
 seqs = numpy.ndarray( int(4e4), dtype='U15000' ) #15 kb
 count = 0

 data = numpy.genfromtxt(filename, dtype='str', delimiter='\n')

 for i in range(0, len(data)):
  if '>' in data[i]:
   names[count] = data[i]
   count = count + 1
  else:
   seqs[count-1] = seqs[count-1] + data[i]

 names = names[0:count]
 seqs = seqs[0:count]

 f = open( filename[:filename.index('.fa')] + '_unwrap' + '.fa', 'w' )
 for i in range(0,len(names)):
  a = names[i]
  b = seqs[i]
  temp = f.write( a + "\n" + b + "\n" )
 
 return f.close()
