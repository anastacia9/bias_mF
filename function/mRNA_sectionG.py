# input:  [0] specify function to be used for calculating the deltaG of a section of interest; options include:
#             [a] "sec_frac"    the deltaG of each partial structure within the section of interest is computed as:
#                               deltaG of whole structure * (# structure nucleotide intervals within section of interest / # nucleotide intervals composing section of interest)
#                               sectional deltaG is the sum of the above deltaG values as well as those of full structures in the section
#             [b] "struc_frac"  the deltaG of each partial structure within the section of interest is computed as:
#                               deltaG of whole structure * (# structure nucleotide intervals within section of interest / # nucleotide intervals composing structure)
#                               sectional deltaG is the sum of the above deltaG values as well as those of full structures in the section
#             [c] "riboG"       an attempt to measure the amount of work a ribosome must do to linearize the mRNA section
#                               structures with positive deltaG require no work to pry apart, so their deltaG is ignored
#                               consider that the ribosome footprint is 30 nucleotides long and structures pried apart once may need to be pried apart a second time
#                               ribosome does not need to unwind the 3'UTR specifically; termination of translation is at the stop codon
#         [1] .csv file with gene names (1st column) and lengths of 5'UTR (2nd column), CDS (3rd column), and 3'UTR (4th column)
#         [2] .csv file with domain location predictions for each gene; first column is gene names (a gene with several domains will appear several times); nucleotide position of domain start in CDS (second column); nucleotide position of domain end in CDS (third column)
#         [3] stdout output .txt file from RNAfold
#         [4] stdout output .txt file from RNAeval
# output: [0] array of gene names
#         [1] array of deltaG values for first 10 nucleotides of 5'UTR
#         [2] array of deltaG values for region -9 to +3 of A[+1] in start codon
#         [3] array of deltaG values for region +4 to +10 of A[+1] in start codon
#         [4] array of deltaG values for stop codon and first 15 nucleotides of 3'UTR
#         [5] array of deltaG values for concatenated sequence remaining when the above four regions are excluded
#         [6] array of deltaG values for concatenated sequence coding for protein domain(s) plus the coding region trailing last domain
#         [7] array of deltaG values for concatenated sequence coding for protein-domain-linker(s) plus the coding region preceding first domain
#         [8] array of deltaG values for whole gene sequence (the 5'UTR + CDS + 3'UTR)
#         [9] array of deltaG values for 5'UTR
#        [10] array of deltaG values for CDS
#        [11] array of deltaG values for 3'UTR
#
# note: an output deltaG value of 1 indicates an absence of sectional deltaG information for that gene
#

def mRNA_sectionG(calcG_function, length_csvfile, domain_csvfile, RNAfold_outputfile, RNAeval_outputfile):

 if "numpy" not in dir():
  import numpy

 if "pandas" not in dir():
  import pandas

 #### Define function for computing deltaG from an array of keys.
 if calcG_function == "sec_frac":            #
  def getG(spec_deltaGs):                    #
   uniqueG = numpy.unique(spec_deltaGs, return_counts=True)
   num_spec_deltaGs = len(spec_deltaGs)      # how many nucleotide intervals long is our section of interest?
   values = 0                                #
   for value in range(0, len(uniqueG[0])):   # for each unique spec_deltaGs key . . 
    key  = uniqueG[0][value]                 # extract unique key
    freq = uniqueG[1][value]                 # extract frequency of unique key
    in_G = key.rfind('_G')                   # find index of deltaG value
    in_s = key.rfind('_s')                   # find index of structure-size value
    G    = float(key[(in_G+2):in_s])         # extract deltaG value
    size = float(key[(in_s+5): ])            # extract structure-size value
    if size == freq:                         # if so, the full structure is contained in our section of interest
     values += G                             # add full deltaG value to 'values' counter
    else:                                    # if so, we have a partial structure contained in our section of interest
     values += G * freq/num_spec_deltaGs     # 
   return values                             # 

 elif calcG_function == "struc_frac":        #
  def getG(spec_deltaGs):                    #
   uniqueG = numpy.unique(spec_deltaGs, return_counts=True)
   num_spec_deltaGs = len(spec_deltaGs)      # how many nucleotide intervals long is our section of interest?
   values = 0                                #
   for value in range(0, len(uniqueG[0])):   # for each unique spec_deltaGs key . . 
    key  = uniqueG[0][value]                 # extract unique key
    freq = uniqueG[1][value]                 # extract frequency of unique key
    in_G = key.rfind('_G')                   # find index of deltaG value
    in_s = key.rfind('_s')                   # find index of structure-size value
    G    = float(key[(in_G+2):in_s])         # extract deltaG value
    size = float(key[(in_s+5): ])            # extract structure-size value
    values += G * freq/size                  # 
   return values                             # 

 elif calcG_function == "riboG":             #
  def getG(spec_deltaGs):                    #
   values = 0                                # tally of negative deltaG values
   for g in range(0, len(spec_deltaGs)):     # for each spec_deltaGs key . . 
    key = spec_deltaGs[g]                    #
    count = g if g < 29 else 29              #
    if key not in spec_deltaGs[g-count:g]:   # if the structure was not just broken by the ribosome
     in_G = key.rfind('_G')                  # find index of deltaG value
     in_s = key.rfind('_s')                  # find index of structure-size value
     G    = float(key[(in_G+2):in_s])        # extract deltaG value
     if G < 0:                               # if deltaG is less than 0 (meaning that the ribosome has to do work to pry the structure apart)
      values += G                            # add deltaG to values counter
   return values                             #
 
 else:
  print('calcG_function not recognized')

 #### Read in stdout output file from RNAfold. Define array of gene names.
 RNAfold     = numpy.genfromtxt(RNAfold_outputfile, dtype='str', delimiter='\n')
 Y_index     = RNAfold[0].find('Y')
 gene_names  = numpy.array([ gene[Y_index:] for gene in RNAfold[0::6] ])
 num_genes   = len(gene_names)

 #### Read in stdout output file from RNAeval.
 RNAeval     = open(RNAeval_outputfile).readlines()

 #### Compute array of gene indices in RNAeval output file (since gene info is not labelled by gene name in this file).
 index_eval  = numpy.array([ line for line in range(0, len(RNAeval)) if '.' in RNAeval[line] ]) + 1
 index_eval  = numpy.append([0], index_eval)

 #### Read in length_csvfile into pandas dataframe.
 df_length   = pandas.read_csv(length_csvfile, index_col = 0)

 #### Read in domain_csvfile into pandas dataframe. Define array of genes with both RNAeval and domain predictions.
 df_domain   = pandas.read_csv(domain_csvfile, index_col = 0)
 temp_index  = numpy.in1d(gene_names, numpy.unique(df_domain.index))
 dgene_names = gene_names[temp_index]

 #### Initialize arrays of gene and gene-region-specific deltaG values. One element per gene.
 gene_deltaG = numpy.ones(num_genes, dtype='float') * 1e6
 cap_deltaG  = numpy.ones(num_genes, dtype='float') * 1e6
 upA_deltaG  = numpy.ones(num_genes, dtype='float') * 1e6
 dwA_deltaG  = numpy.ones(num_genes, dtype='float') * 1e6
 dwS_deltaG  = numpy.ones(num_genes, dtype='float') * 1e6
 EvE_deltaG  = numpy.ones(num_genes, dtype='float') * 1e6
 Do3_deltaG  = numpy.ones(num_genes, dtype='float') * 1e6
 Li5_deltaG  = numpy.ones(num_genes, dtype='float') * 1e6
 UT5_deltaG  = numpy.ones(num_genes, dtype='float') * 1e6
 CDS_deltaG  = numpy.ones(num_genes, dtype='float') * 1e6
 UT3_deltaG  = numpy.ones(num_genes, dtype='float') * 1e6

 #### For each gene . . .                                              #
 for i in range(0, num_genes):                                         #
     info       = RNAeval[index_eval[i] : index_eval[i+1]-2]           # extract structure predictions
     lengths    = df_length.loc[gene_names[i]]                         # extract gene 5'UTR, CDS, and 3'UTR lengths
     len5       = lengths['UTR5_len']                                  # define 5'UTR length
     lenc       = lengths['CDS_len']                                   # define CDS length
     len3       = lengths['UTR3_len']                                  # define 3'UTR length
     length     = sum(lengths)                                         # compute 5'UTR + CDS + 3'UTR length
     deltaGs    = numpy.ndarray(length, dtype='U150')                  # initialize empty deltaG array (one element per nucleotide interval plus an extra element at index=0)
     multi_size = numpy.zeros(length, dtype='float')                   # initialize multi-loop size array
     multi_loc  = numpy.array([None]*length)                           # initialize multi-loop location array
     multi_key  = numpy.ndarray(length, dtype='U150')                  # initialize multi-loop deltaG-key array
     num_multi  = 0                                                    # counter for multi-loops
     G_ext_key  = None                                                 # marks presence of absence of external loop
                                                                       #
     ### For the each of the structures in the mRNA . . .              #
     for structure in info:                                            #
       index = structure.rfind( ' ' ) + 1                              # find index of deltaG value in structure information
       deltaG = float(structure[index:])                               # extract the deltaG value
       a0     = structure.find( ')' )                                  # find the first base-pair location of structure                                                              #
       a1     = structure[15:a0].split( ',' )                          #
       if 'Interior' in structure:                                     #
         paira = list(map(int, a1))                                    # save first base-pair location
         b0    = structure.rfind( '(' ) + 1                            # find the second base-pair location of interior loop
         b1    = structure.rfind( ')' )                                #
         b2    = structure[b0:b1].split(',')                           #
         pairb = list(map(int, b2))                                    # save second base-pair location
         diff1 = pairb[0] - paira[0]                                   # compute distance between nucleotides on one side of loop
         diff2 = paira[1] - pairb[1]                                   # compute distance between nucleotides on other side of loop
         size  = diff1 + diff2                                         # compute total size of interior loop
         key   = 'interior_'+str(paira[0])+'_'+str(paira[1])+'_G'+str(deltaG)+'_size'+str(size)
         deltaGs[paira[0]:pairb[0]] = key                              # add interior loop deltaG key to deltaG array
         deltaGs[pairb[1]:paira[1]] = key                              #
       elif 'Hairpin' in structure:                                    #
         paira = list(map(int, a1))                                    #
         size  = paira[1] - paira[0]                                   # compute size of hairpin loop
         deltaGs[paira[0]:paira[1]] = 'hairpin_'+str(paira[0])+'_'+str(paira[1])+'_G'+str(deltaG)+'_size'+str(size)
       elif 'Multi' in structure:                                      #
         paira = list(map(int, a1))                                    #
         multi_size[num_multi] = paira[1]-paira[0]                     # store size of multi loop
         multi_loc[num_multi]  = paira                                 # store range of multi loop
         multi_key[num_multi]  = 'multi_'+str(paira[0])+'_'+str(paira[1])+'_G'+str(deltaG)
         num_multi += 1                                                # count multi loop
       elif 'External' in structure:                                   #
         G_ext_key = '_G' + str(deltaG)                                #
                                                                       #
     ### Work out multi loops of mRNA.                                 #
     multi_size  = multi_size[:num_multi]                              #
     multi_index = numpy.argsort(multi_size)                           # sort multi loops from smallest to largest and store indices
     multi_loc   = multi_loc[:num_multi][multi_index]                  #
     multi_key   = multi_key[:num_multi][multi_index]                  #
     for j in range(0, num_multi):                                     # for each multi loop . . .
      fraction = deltaGs[multi_loc[j][0]:multi_loc[j][1]]              # extract elements corresponding to multi-loop range in deltaGs array
      where_0  = numpy.where(fraction == '')[0]                        # which of these elements is not yet defined with a deltaG key?
      fraction[where_0] = multi_key[j]+'_size'+str(len(where_0))       # append size of where_0 to multi loop key; populate deltaGs array with key
                                                                       #
     ### Work out external loops of mRNA.                              #
     if G_ext_key != None:                                             # add external loop data to deltaGs array (if applicable)
      where_0 = numpy.where(deltaGs == '')[0]                          # extract empty elements of deltaGs array
      deltaGs[where_0] = G_ext_key+'_size'+str(len(where_0)-1)         # populate empty elements with external loop key
      deltaGs = deltaGs[1:]                                            # get rid of first element in deltaGs array (since RNAeval indexes from 1)
                                                                       #
     ### Define matching array of deltaGs for each section.            #
     drop = 0                                                          #
     gene_deltaG[i] = getG(deltaGs) if calcG_function != "riboG" else getG(deltaGs[:len5+lenc+15-1])
     UT5_deltaG[i]  = getG(deltaGs[:len5-1])                           # define deltaG of 5'UTR
     CDS_deltaG[i]  = getG(deltaGs[len5:len5+lenc-1])                  # define deltaG of CDS
     UT3_deltaG[i]  = getG(deltaGs[len5+lenc:])                        # define deltaG of 3'UTR
     dwA_deltaG[i]  = getG( deltaGs[(len5+4-1):(len5+10-1)] )          # define deltaG of region downstream of start codon
     if len5 >= 19:                                                    # we have a section overlap problem if the 5'UTR is less than 19 nucleotides long
       cap_deltaG[i] = getG( deltaGs[0:9] )                            # define deltaG of cap region
       upA_deltaG[i] = getG( deltaGs[(len5-9):(len5+3-1)] )            # define deltaG of region upstream and including start codon
       drop += 1                                                       #
     if len3 != 0:                                                     # we have a problem if the 3'UTR is of length zero
       dwS_deltaG[i] = getG( deltaGs[(len5+lenc-3):(len5+lenc+15-1)] ) # define deltaG of region downstream and including stop codon
       drop += 1                                                       #
     if drop == 2:                                                     #
      EvE_deltaG[i] = gene_deltaG[i]-cap_deltaG[i]-upA_deltaG[i]-dwA_deltaG[i]-dwS_deltaG[i]
     if gene_names[i] in dgene_names:                                  #
       keep = 0                                                        #
       domain_info = df_domain.loc[gene_names[i]]                      # extract domain location information (relative to CDS)
       vals        = domain_info.get_values().flatten()                # flatten array of domain coordinates
       num_vals    = len(vals)                                         #
       sorted      = numpy.argsort(vals[0::2])                         # sort array of domain start coordinates
       intervals   = numpy.zeros(num_vals, dtype='int')                #
       intervals[0::2] = vals[0::2][sorted]                            #
       intervals[1::2] = vals[1::2][sorted]                            #
       if num_vals == sum(numpy.sort(intervals) == intervals) and num_vals == len(set(intervals)): # check that domain intervals do not overlap AND check that the end of one domain is not at the same position as the start of the following domain
        if lenc-intervals[-1] > 1:                                     # is there a 3' coding region more than 1 nucleotide in length?
         Do3_deltaG[i] = getG(deltaGs[len5+intervals[-1]:len5+lenc-1]) # calculate deltaG of 3' coding sequence
         keep = 1                                                      # specify that deltaG for 3' region was extracted
        for k in range(0, num_vals, 2):                                # for each domain . . .
         Do3_deltaG[i] = 0 if k==0 and keep!=1 else Do3_deltaG[i]      # conditional for 3' coding region deltaG
         first = len5+intervals[k] - 1                                 #
         last  = len5+intervals[k+1] - 1                               #
         Do3_deltaG[i] += getG( deltaGs[first:last] )                  # calculate deltaG of each domain and sum together
        if intervals[0] > 2:                                           # is there a 5' coding region more than 1 nucleotide in length?
         Li5_deltaG[i] = getG( deltaGs[len5:len5+intervals[0]-2] )     # calculate deltaG of 5' coding sequence
         keep = 2                                                      #
        if num_vals != 2:                                              # if there are linker regions
          for k in range(0, num_vals-2, 2):                            # for each linker . . .
           first = len5+intervals[k+1]                                 # 
           last  = len5+intervals[k+2] - 2                             #
           if first < last:                                            # if the linker region is more than 1 nucleotide in length
            Li5_deltaG[i] = 0 if k==0 and keep!=2 else Li5_deltaG[i]   #
            Li5_deltaG[i] += getG(deltaGs[first:last])                 # calculate deltaG of each linker and sum together
 return gene_names, cap_deltaG, upA_deltaG, dwA_deltaG, dwS_deltaG, EvE_deltaG, Do3_deltaG, Li5_deltaG, gene_deltaG, UT5_deltaG, CDS_deltaG, UT3_deltaG
