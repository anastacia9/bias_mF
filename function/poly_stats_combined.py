# How many polymorphisms in combined gene regions (D+3' and L+5')?
# How many unique variants per polymorphism?
# How widespread is the most frequenct variant?
# input:  [0] combined region type - options are 'domain+threeprime' and 'linker+ramp'
#         [1] .csv file with domain predictions (cannot include genes outside of the 1636 we have protein:mRNA data for)
#         [2] path to fasta files with gene CDS sequences; gene order must be identical across fasta files
# output: [0] array of gene names
#         [1] array of combined region lengths (one entry per gene)
#         [2] array of polymorphism counts in combined region (one entry per gene across strains)
#         [3] array of median number of variants per polymorphism in the combined region (one entry per gene across strains)
#         [4] array of mean   number of variants per polymorphism in the combined region (one entry per gene across strains)
#         [5] array of mean maximum variant frequency per polymorphism in the combined region (one entry per gene across strains)
#         [6] array of polymorphism locations in the combined region

def poly_stats_combined(combined_region_type, csv_file, path):

 if "os" not in dir():
   import os

 if "numpy" not in dir():
  import numpy

 if "pandas" not in dir():
  import pandas

 if "dltAI_mod" not in dir():
  import dltAI_mod

 if combined_region_type == 'domain+threeprime':
  x = 7
  y = 10
  drops = ['YDR150W', 'YGR256W', 'YKR103W', 'YBR096W', 'YDR009W', 'YNR017W', 'YHR089C', 'YOL060C', 'YLR345W', 'YNR069C', 'YML072C', 'YPR091C', 'YCR098C', 'YJL104W', 'YML052W', 'YLR303W', 'YGR100W', 'YHR163W', 'YPR041W', 'YDL122W']

 if combined_region_type == 'linker+ramp':
  x = 8
  y = 9
  drops = ['YDR150W', 'YKR103W', 'YBR096W', 'YDR009W', 'YNR069C', 'YML072C', 'YJL104W', 'YHR163W']

 # Define path to strain sequence files.
 files = os.listdir(path)
 num_files = len(files)

 # Run dltAI_mod.py on each file and store results in array.
 dltAI_array = [None] * num_files
 count = 0
 for filename in files:
  dltAI_array[count] = dltAI_mod.dltAI(path + '\\' + filename, csv_file, 'tAI') 
  count += 1
  print(filename)

 # Define vector of gene names.
 gene_names = dltAI_array[0][0]
 num_genes = len(gene_names)

 # Define empty arrays.
 region_length = numpy.zeros(num_genes, dtype='int')
 num_poly = numpy.zeros(num_genes, dtype='int')
 med_poly = numpy.zeros(num_genes, dtype='int')
 mean_poly = numpy.zeros(num_genes, dtype='int')
 prop_max_var = numpy.zeros(num_genes)
 loc_poly = [None] * num_genes


 # Calculate!
 for i in range(0, num_genes):
  if gene_names[i] in drops:
   continue
  region1_seq = dltAI_array[0][x][i]
  region2_seq = dltAI_array[0][y][i]
  if region1_seq or region2_seq:
   region_length[i] = len( ''.join(region1_seq)+region2_seq )
   df = pandas.DataFrame(columns=numpy.arange(region_length[i]), index=numpy.arange(num_files))
   for j in range(0, num_files):
    region_seq = ''.join(dltAI_array[j][x][i]) + dltAI_array[j][y][i]
    df.loc[j] = [region_seq[k] for k in range(0, region_length[i])]
   num_variants_total = numpy.array(list( df.nunique(axis=0)-1 ))
   loc_poly[i] = numpy.where(num_variants_total > 0)[0]
   num_poly[i] = len(loc_poly[i])
   if num_poly[i]:
    med_poly[i] = numpy.median(num_variants_total[loc_poly[i]])
    mean_poly[i] = numpy.mean(num_variants_total[loc_poly[i]])
    prop_max_var[i] = numpy.mean( numpy.array([max(df[m].value_counts(normalize=True)) for m in loc_poly[i]]) )

 return gene_names, region_length, num_poly, med_poly, mean_poly, prop_max_var, loc_poly, dltAI_array

