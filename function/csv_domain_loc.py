# input:  [0] domain assignments file for individual genes
#             downloaded from http://downloads.yeastgenome.org/curation/calculated_protein_info/domains/domains.tab
#             and based on InterProScan on January 14, 2014
#         [1] domain database (e.g. 'Pfam', 'PANTHER', 'Gene3D', 'PRINTS', ...)
# output: [0] 'domains.csv' file with domain location predictions for each gene as specified by the input domain database
#

def csv_domain_loc(assignments_file, domain_database):

 if "numpy" not in dir():
  import numpy

 if "pandas" not in dir():
  import pandas

 data = numpy.genfromtxt(assignments_file, dtype='str', delimiter='\n')
 num_entries = len(data)
 gene_names = numpy.ndarray(num_entries, dtype='U15')
 domain_starts = numpy.zeros(num_entries, dtype=int)
 domain_ends = numpy.zeros(num_entries, dtype=int)

 count = 0
 for i in data:
  entry = i.split('\t')
  if domain_database in entry:
   gene_names[count] = entry[0]
   domain_starts[count] = int(entry[6])
   domain_ends[count] = int(entry[7])
   count += 1

 gene_names = gene_names[:count]
 domain_starts = domain_starts[:count]
 domain_ends = domain_ends[:count]

 index = numpy.argsort(gene_names)
 gene_names = gene_names[index]
 domain_starts = domain_starts[index]
 domain_ends = domain_ends[index]

 df = pandas.DataFrame(numpy.array([domain_starts, domain_ends]).T, index=gene_names)
 df.columns = ['domain_start', 'domain_end']

 df.to_csv('domains.csv', )
   
