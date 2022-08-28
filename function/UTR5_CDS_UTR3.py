# input:  [0] .fa filename with chromosome DNA sequences (filename in single quotes); filename wrapped at 60 characters
#         [1] .fa filename with gene DNA sequences (filename in single quotes); filename wrapped at 60 characters
#         [2] .csv filename with [a] gene name in first column
#                                [b] chromosome number in second column (e.g. 'chrI', 'chrII', 'chrIII', ...)
#                                [c & d] leading strand chromosome coordinates (position1 and position2) of CDS in third and fourth columns, respectively
#                                [e] leading (+) or lagging (-) strand specification in fifth column
#                                [f] 5'UTR length in sixth column
#                                [g] 3'UTR length in seventh column
# output: [0] array of gene names
#         [1] array of corresponding 5'UTR sequences (one entry per gene)
#         [2] array of corresponding CDS sequences (one entry per gene)
#         [3] array of corresponding 3'UTR sequences (one entry per gene)
#         [4] array of corresponding concatenated 5'UTR, CDS, and 3'UTR (in that order) sequences (one entry per gene)

def UTR5_CDS_UTR3(chr_filename, gene_filename, csv_filename):
 
 if "numpy" not in dir():
  import numpy

 if "pandas" not in dir():
  import pandas

 def rev_transcribe(input_sequence):
  replace = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
  seq_array = numpy.array([ replace[letter] for letter in input_sequence ])
  return ''.join(seq_array)

 data_chr   = numpy.genfromtxt(chr_filename, dtype='str', delimiter='\n')
 data_gene  = numpy.genfromtxt(gene_filename, dtype='str', delimiter='\n')
 df_UTR_loc = pandas.read_csv(csv_filename, index_col=0)

 gene_names   = numpy.array(list( df_UTR_loc.index ))
 num_genes    = len(gene_names)
 UTR5_array   = numpy.ndarray(num_genes, dtype='U30000')
 CDS_array    = numpy.ndarray(num_genes, dtype='U30000')
 UTR3_array   = numpy.ndarray(num_genes, dtype='U30000')
 concat_array = numpy.ndarray(num_genes, dtype='U30000')

 chr_indices = numpy.array([0, 3838, 17393, 22671, 48205, 57821, 62325, 80509, 89888, 97221, 109652, 120767, 138738, 154147, 167221, 185411, 201214])
 chr_names   = numpy.array(['chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII', 'chrIX', 'chrX', 'chrXI', 'chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI', 'chrM'])
 chr_dict    = dict(zip(chr_names, chr_indices))

 for i in range(0, num_genes):
  info = df_UTR_loc.iloc[i]
  where_chr = chr_dict[info['chr']] # index of chromosome in data_chr
  if info['lead_or_lag'] == '+':    # information is on leading or lagging strand?
   A = 'UTR5'
   B = 'UTR3'
  else:
   A = 'UTR3'
   B = 'UTR5'
  lenA = info[A]
  lenB = info[B]
  lengths = numpy.array([lenA, lenA, lenB, lenB])
  UTRA_chr_pos1  = info['position1'] - info[A]  # relative to the chromosome, where on the leading strand does the A UTR start?
  UTRA_chr_pos2  = info['position1'] - 1        # relative to the chromosome, where on the leading strand does the A UTR end?
  CDS_chr_pos1   = info['position1']            # relative to the chromosome, where on the leading strand does this gene start?
  CDS_chr_pos2   = info['position2']            # relative to the chromosome, where on the leading strand does this gene end?
  UTRB_chr_pos1  = info['position2'] + 1        # relative to the chromosome, where on the leading strand does the B UTR start?
  UTRB_chr_pos2  = info['position2'] + info[B]  # relative to the chromosome, where on the leading strand does the B UTR end?
  region_indices = numpy.array([UTRA_chr_pos1, UTRA_chr_pos2, UTRB_chr_pos1, UTRB_chr_pos2])
###
  CDS_index = numpy.where(data_gene == '>'+gene_names[i])[0][0] + 1
  CDS_array[i] = data_gene[CDS_index]
###
  for region in numpy.array([0, 2]):
   ###initial section###
   pos1 = region_indices[region]                        # relative to the chromosome, where on the leading strand does this region start?
   pos2 = region_indices[region+1]                      # relative to the chromosome, where on the leading strand does this region end?
   length = lengths[region]                             # how long is the region?
   groups60_to_pos1 = int(numpy.floor( (pos1-1)/60 ))   # how many groups of 60 nucleotides precede pos1?
   left_over_to_pos1 = (pos1-1) % 60                    # after how many nucleotides into the following group of 60 is pos1?
   pos1_post = data_chr[where_chr + groups60_to_pos1 + 1][(left_over_to_pos1):]
   ##
   if len(pos1_post) < length:                          # is pos1_post length less than the region length?
    ###final section###
    groups60_to_pos2 = int(numpy.floor( (pos2-1)/60 ))  # how many groups of 60 nucleotides precede pos2?
    left_over_after_pos2 = (pos2-1) % 60                # after how many nucleotides into the following group of 60 is pos2?
    pos2_pre = data_chr[where_chr + groups60_to_pos2 + 1][:(left_over_after_pos2+1)]
    ###middle section###
    middle = data_chr[(where_chr+groups60_to_pos1+2):(where_chr+groups60_to_pos2+1)]
    middle_join = ''.join(middle)
   ##
   else:
    pos1_post = pos1_post[:length]
    middle_join = ''
    pos2_pre = ''
   ##
   ###populate arrays###
   if region == 0:
    if info['lead_or_lag'] == '-':
     temp0 = ''.join([pos1_post, middle_join, pos2_pre])
     temp1 = temp0[::-1]
     UTR3_array[i] = rev_transcribe(temp1)
    else:
     UTR5_array[i] = ''.join([pos1_post, middle_join, pos2_pre])
   if region == 2:
    if info['lead_or_lag'] == '-':
     temp2 = ''.join([pos1_post, middle_join, pos2_pre])
     temp3 = temp2[::-1]
     UTR5_array[i] = rev_transcribe(temp3)
    else:
     UTR3_array[i] = ''.join([pos1_post, middle_join, pos2_pre])
    concat_array[i] = ''.join([UTR5_array[i], CDS_array[i], UTR3_array[i]])


 return gene_names, UTR5_array, CDS_array, UTR3_array, concat_array

