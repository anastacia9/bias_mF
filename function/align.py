# input: [0] path to folder name containing unwrapped and translated fasta files with gene names organized alphabetically (path in single quotes); no limit to number of files; files must be present in working directory
#        [1] number of total genes to be aligned
# output: wrapped file "alignments.fa" in working directory; progress reports printed for every 100 genes
#
# note: MUSCLE provided by Edgar, R.C. Nucleic Acids Res 32(5), 1792-97;
#       strains must have the same set of genes (the 22 strains from YRC do);

def align(path, num_genes):

 if "os" not in dir():
  import os

 if "numpy" not in dir():
  import numpy

 directory = os.fsencode(path)

 print('number of genes aligned:')
 count = 0
 for i in range(0, num_genes):
  h = open( 'temp.fa', 'w' )
  for file in os.listdir(directory):
   filename = os.fsdecode(file)
   if filename.endswith(".fa"):
    with open(filename) as f:
     lines = f.readlines()
    line = lines[2*i]
#    temp = h.write(  line + lines[2*i+1] )
    temp = h.write(  line[0:line.find('/n')] + '   ' + filename + '\n' + lines[2*i+1] )
  h.close()
  temp1 = os.system('muscle.exe -in temp.fa >>alignments.afa 2>nul')
  count = count + 1
  if count % 100 == 0:
   print(count)
 temp2 = os.remove('temp.fa')


 return print('complete')
