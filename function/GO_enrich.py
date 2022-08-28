# input:  [0] list of file names (one file per set of genes in question)
#             file can be generated via the GO Enrichment Analysis webtool at geneontology.org
#             each file should be tab-delimited and contain the following columns (with header row)
#             [a] GO biological process term
#             [b] number of genes in the full yeast genome that are associated with the corresponding GO term
#             [c] number of genes in the subset of interest that are associated with the corresponding GO term
#             [d] number of genes in the subset of interest expected to be associated with the corresponding GO term
#             [e] fold enrichment of corresponding GO term in subset of interest ('< 0.01' values will be converted to 10e-5)
#             [f] +/- enrichment specification
#             [g] raw p-value
#             [h] FDR-value
#         [1] FDR-value cutoff as float (i.e. ignore GO terms with FDR-value greater than this cutoff)
#             if no cutoff, enter "1"
# output: [0] array of GO-term arrays; one array (per file in input list) of GO-terms ordered from low to high FDR-value
#         [1] array of gene count arrays; one array (per file in input list) of gene counts ordered from low to high FDR-value
#         [2] array of fold enrichment value arrays; one array (per file in input list) of fold enrichment values ordered from low to high FDR-value
#         [3] array of FDR-value arrays; one array (per file in input list) of FDR-values ordered from low to high FDR-value
#
# my_list = ['GO_enrichment_1620.txt', 'GO_enrichment_1458.txt', 'GO_enrichment_1092.txt', 'GO_enrichment_779.txt', 'GO_enrichment_lowtAI723.txt', 'GO_enrichment_hightAI724.txt', 'GO_enrichment_lowG723.txt', 'GO_enrichment_highG724.txt', 'GO_enrichment_185.txt']

def GO_enrich(my_list, val):

 if "numpy" not in dir():
  import numpy
 
 num_files = len(my_list)
 data = [numpy.genfromtxt(file, dtype='str', delimiter='\n') for file in my_list]
 terms = [None] * num_files
 genes = [None] * num_files
 folds = [None] * num_files
 FDRs  = [None] * num_files

 count = 0
 for datum in data:
  num_terms = len(datum)
  temp_terms = numpy.ndarray(num_terms, dtype='U100')
  temp_genes = numpy.zeros(num_terms, dtype='int')
  temp_folds = numpy.zeros(num_terms, dtype='float')
  temp_FDRs  = numpy.zeros(num_terms, dtype='float')
  for i in range(1, num_terms):
   info = datum[i].split('\t')
   temp_terms[i] = info[0]
   temp_genes[i] = int(info[2])
   if info[4] == '< 0.01 ':
    info[4] = 10e-5
   temp_folds[i] = float(info[4])
   temp_FDRs[i]  = float(info[7])
  loc = numpy.where(temp_FDRs[1:] < val)[0] # where are FDR values less than val?
  srt = numpy.argsort(temp_FDRs[1:][loc])   # sort (by index) the FDR values defined by loc in descending order
  terms[count] = temp_terms[1:][loc][srt]
  genes[count] = temp_genes[1:][loc][srt]
  folds[count] = temp_folds[1:][loc][srt]
  FDRs[count]  = temp_FDRs[1:][loc][srt]
  count += 1

 return terms, genes, folds, FDRs



