# input: [0] fasta file with gene sequence information for one strain; (filename in single quotes); genes must be in alphabetical order;
#        [1] .csv file with transcript abundance data; first column has gene names and successive column(s) labelled with strain name(s);
#            genes must be in alphabetical order; (csv filename in single quotes);
#        [2] name of strain in question;
#        [3] .gff general feature format file with gene descriptions including sequence type; (filename in single quotes);
#
# output: [0] dictionary of nTE values (one value for each codon) to be used for ntAI calculations;
#
# Note: codoncount.py and codontable_tAI.py must be in working directory

def codontable_ntAI(filename, transcripts_csv, strain_name, gff_file):

 if "codoncount" not in dir():
  import codoncount

 if "tAI_cTE_codontable.py" not in dir():
  import tAI_cTE_codontable

 if "pandas" not in dir():
  import pandas

 if "numpy" not in dir():
  import numpy

 if "gffgenenames" not in dir():
  import gffgenenames

# Load transcripts .csv file into pandas dataframe.
 dft = pandas.read_csv(transcripts_csv, index_col=0)
 transcript_data = numpy.array(list(dft[strain_name]))

# Define gene names.
 gene_names = numpy.array(list(dft.index))
 gene_names = numpy.array(['>' + name for name in gene_names])
 gene_counts = len(gene_names)

# Load fasta file with gene sequence information for a single strain.
 sequences = numpy.genfromtxt(filename, dtype='str', delimiter='\n')
 seq_titles = sequences[0::2]

# Count codon frequencies in each gene with codoncount.py.
 freqs = codoncount.codoncount(filename)
 index = numpy.where( numpy.in1d(seq_titles, gene_names) )[0]
 if sum(seq_titles[index] == gene_names) != gene_counts:
  print('invalid result; make sure that genes are in alphabetical order in fasta file and in csv file')
 counts = freqs[0][index]

# Compute codon usage.
 cu = sum(  numpy.array([ transcript_data[i]*counts[i] for i in range(0, gene_counts) ])  )
 cu = cu/max(cu)
 dict_codons = dict(zip(freqs[2], cu))

# Compute cTE and then nTE values.
 cTE = tAI_cTE_codontable.tAI_cTE_codontable(gff_file)
 temp = cTE[0]
 dict_w = {k: temp[k]/dict_codons[k] for k in dict_codons}
 my_keys = list( dict_w.keys() )
 my_values = list( dict_w.values() )
 my_values /= max(my_values)
 dict_w = dict(zip(my_keys, my_values))

 return dict_w

