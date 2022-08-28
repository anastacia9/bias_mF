# input:  [0] file with coding sequences in FASTA format (file name in single quotes);
#         [1] .csv file with domain location predictions for each gene (output of 'csv_domain_loc.py');
#            first column is gene name;
#            second column is nucleotide position of domain start;
#            third column is nucleotide position of domain end;
#         [2] codon relative adaptiveness table in dictionary format; if default (dos Reis et al, 2003), type 'tAI';
#
# output: [0] array of gene names
#
#         [1] array of tAI values for gene-level concatenated domain sequences; one entry per gene
#         [2] array of tAI values for gene-level concatenated linker sequences; one entry per gene
#         [3] array of rare codon ramp (coding section preceding first domain) tAI values; one entry per gene
#         [4] array of 3' (coding section trailing last domain) tAI values; one entry per gene
#
#         [5] array of arrays with tAI calculated for individual gene domain sequences; one array per gene
#         [6] array of arrays with tAI calculated for individual gene linker sequences; one array per gene
#
#         [7] array of arrays with individual gene domain sequences; one array per gene
#         [8] array of arrays with individual gene linker sequences; one array per gene
#         [9] array of rare codon ramp sequences; one entry per gene
#        [10] array of 3' sequences; one entry per gene

def dltAI(sequence_file, csv_file, codontable):

 if "tai_seq" not in dir():
  import tai_seq

 if "numpy" not in dir():
  import numpy

 if "pandas" not in dir():
  import pandas

###
# Define input sequence file as numpy array.
 sequences = numpy.genfromtxt(sequence_file, dtype='str', delimiter='\n')

###
# Define input .csv file with domain locations as pandas dataframe.
 df = pandas.read_csv(csv_file)

###
# Create an array of domain intervals for each gene.
 all_genes = numpy.array( list(df['gene']) )
 genes = numpy.unique(all_genes)
 all_starts = numpy.array( list(df['domain_start']) )
 all_ends = numpy.array( list(df['domain_end']) )

 num_genes = len(genes)
 interval_values = [None] * num_genes

###
# Define empty arrays.
 cdtAI =           [None] * num_genes
 cltAI =           [None] * num_genes
 ramptAI =         [None] * num_genes
 three_primetAI =  [None] * num_genes
 dtAI =            [None] * num_genes
 ltAI =            [None] * num_genes
 dseqs =           [None] * num_genes
 lseqs =           [None] * num_genes
 rampseqs =        [None] * num_genes
 three_primeseqs = [None] * num_genes
 cdseqs =          [None] * num_genes
 clseqs =          [None] * num_genes

###
### Perform tAI computations, one gene at a time.
###

###
# Find index of gene in sequences file.
 for count in range(0, num_genes):
  gene_index = numpy.where(sequences == '>' + genes[count])[0][0]

###
# Create array of domain interval values.
  index = numpy.where(all_genes == genes[count])[0]
  num_intervals = len(index)
  temp1 = all_starts[index]
  temp2 = all_ends[index]
  index = numpy.argsort(temp1)
  interval = numpy.zeros(num_intervals*2, dtype=int)
  interval[0::2] = temp1[index]
  interval[1::2] = temp2[index]
  interval_values[count] = interval

  ##Check that domain intervals do not overlap.
  if sum( numpy.sort(interval_values[count]) == interval_values[count] ) != 2*num_intervals:
   continue

  ##Check that the end of one domain is not at the same position as the start of the following domain.
  if len(interval_values[count]) > len(set(interval_values[count])):
   continue

  interval_values[count] += num_intervals*[-1,0]

###
# Specify sequences. Search for premature stop codons first.
  gene_seq = sequences[gene_index+1]
  codons = [ gene_seq[i:i+3] for i in range(0, len(gene_seq), 3) ]
  stops = ['TGA', 'TAG', 'TAA']
  firststop = numpy.where(numpy.in1d(codons, stops) == True)[0]
  if len(firststop):
   codons = codons[:firststop[0]+1]
  sequence_sections = numpy.split(codons, interval_values[count])
  for i in range(0, len(sequence_sections)):
   sequence_sections[i] = ''.join(sequence_sections[i])

###
# Define domain and concatenated domain sequences and tAI values.
  temp0 = sequence_sections[1::2]
  dseqs[count] = temp0
  cdseqs[count] = ''.join(temp0)
  dtAI[count] = numpy.array([ tai_seq.tAI_seq(j, codontable) for j in temp0 ])
  cdtAI[count] = tai_seq.tAI_seq(cdseqs[count], codontable)

###
# Define linker and concatenated linker sequences and tAI values.
  temp1 = sequence_sections[2:-1:2]
  lseqs[count] = temp1
  clseqs[count] = ''.join(temp1)
  ltAI[count] = numpy.array([ tai_seq.tAI_seq(k, codontable) for k in temp1 ])
  cltAI[count] = tai_seq.tAI_seq(clseqs[count], codontable)

###
# Define ramp sequence and tAI value.
  rampseqs[count] = sequence_sections[0]
  ramptAI[count] = tai_seq.tAI_seq(rampseqs[count], codontable)

###
# Define 3' sequence and tAI value.
  temp2 = sequence_sections[-1]
  if ( len(temp0) + len(temp1) + 1 ) != len(sequence_sections): #if there are 3' sections
   three_primeseqs[count] = temp2
   three_primetAI[count] = tai_seq.tAI_seq(temp2, codontable)

 return genes, cdtAI[:count+1], cltAI[:count+1], ramptAI[:count+1], three_primetAI[:count+1], dtAI[:count+1], ltAI[:count+1], dseqs[:count+1], lseqs[:count+1], rampseqs[:count+1], three_primeseqs[:count+1]
