# This python script provides several options for converting a cell x gene matrix into an h5ad file appropriate for MapMyCells.  It is divided into two parts:
# Part 1: DATA INPUT. This section describes how to read your data in from multiple starting file formats and store as a matrix in python, with associated vectors for gene names and cell names.
# Part 2: DATA OUTPUT. This section describes how to output your variable in a compressed h5ad file, check the size, and then split the output file into multiple files for upload to MapMyCells if the size exceeds 2GB.

#######################################################

# PART 1: DATA INPUT

# The first part describes how to read in your data and store as a matrix. Note that only one section in this part needs to be run!
# Option 1. If your data is stored as a csv file with CELL names as columns and GENE names in the first row (e.g., matrix requires transposition)
# Option 2. If your data is stored as a csv file with GENE names as columns and CELL names in the first row
# Option 3. If your data is stored as an hdf5 file with CELLS as columns and GENES as rows
# Option 4. If your data is already in h5ad format, but you want to compress it.

# At the end of this section you will have a matrix called "count_matrix" a sample vector called "obs" and a gene vector called "var"


### OPTION 1: csv file stored with cells as columns and genes as rows

# As a specific example, we'll download SmartSeq data from mouse VISp from here (https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-v1-and-alm-smart-seq) and convert to h5ad. Prior to running any code, download the data from this file (https://celltypes.brain-map.org/api/v2/well_known_file_download/694413985)[292MB] to your working directory, and unzip it. For simplicity, we will only consider the exons.

# Import libraries
import pandas as pd

# Read in count matrix numbers and store as dataframe
exons = pd.read_csv("mouse_VISp_2018-06-14_exon-matrix.csv") 
exons.rename({exons.columns[0]: 'gene_entrez_id'}, axis=1, inplace=True)
count_matrix = exons[exons.columns[1:]].to_numpy()

# Store sample-based metadata in a dataframe called obs
obs = pd.DataFrame(
     [{'sample_id': s} for s in exons.columns[1:]]).set_index('sample_id')

# Store gene-based metadata (entrez IDs and gene symbols) in a dataframe called var
var = exons[['gene_entrez_id']]
gene_symbol_df = pd.read_csv('mouse_VISp_2018-06-14_genes-rows.csv')
var = var.join(gene_symbol_df.set_index('gene_entrez_id'), on='gene_entrez_id')
var = var.set_index('gene_symbol')

# Transpose counts so samples are rows and genes are columns. Note that we've downsampled in this example since the matrix is huge.
#count_matrix = count_matrix.transpose()  # CODE YOU WOULD ACTUALLY USE
count_matrix = count_matrix[:, :10].transpose()
obs = obs[:10]


### OPTION 2: csv file stored with genes as columns and cells as rows

# For this option, we do not have an example file, but replace "MATRIX.csv" with the actual name of your matrix, and "sample_id" with the actual first column name of your data matrix (e.g., the column name for the cell ids).

# Import libraries
import pandas as pd

# Read in count matrix numbers and store as dataframe
exons = pd.read_csv('MATRIX.csv') 

exons.rename({exons.columns[0]: 'sample_id'}, axis=1, inplace=True) 
count_matrix = exons[exons.columns[1:]].to_numpy()

# Store sample-based metadata in a dataframe called obs
sample_ids = exons['sample_id']

obs = pd.DataFrame(
     [{'sample_id': s} for s in sample_ids]).set_index('sample_id')

# Store gene-based metadata in a dataframe called var
var = pd.DataFrame(
    [{'gene_symbol': g} for g in exons.columns[1:]]).set_index('gene_symbol')
		

### OPTION 3: hdf5 file with cells as columns and genes as rows

# As a specific example, we'll download 10x data from mouse cortex and hippocampus from here (https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x). Prior to running any code, download the data from this file (https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_10x/expression_matrix.hdf5)[5.3GB] to your working directory

# Import libraries
import h5py
import pandas as pd

# Read in count matrix numbers and store as dataframe. Note that we are downsampling before we transpose because this file is huge.
src = h5py.File('expression_matrix.hdf5', 'r')
#count_matrix = src['data/counts'][()].transpose()   # CODE YOU WOULD ACTUALLY USE
count_matrix = src['data/counts'][:, :10].transpose()
sample_names = src['data/samples'][()][:10] # REMOVE [:10] for actual code
gene_symbols = src['data/gene'][()]
src.close()

# Store sample-based metadata in a dataframe called obs
obs = pd.DataFrame(
    {'sample_name': s} for s in sample_names).set_index('sample_name')

# Store gene-based metadata in a dataframe called var
var = pd.DataFrame(
    {'gene_symbol': g} for g in gene_symbols).set_index('gene_symbol')
		

### OPTION 4: If your data is already in h5ad format, but you want to compress it

# No specific example provided here, but replace 'path/to/file/my_data.h5ad' with the actual path to your file.

# Import libraries
import anndata

# Read in the data.
my_adata = anndata.read_h5ad('path/to/file/my_data.h5ad')
count_matrix = my_adata.X

# Store sample-based metadata in a dataframe called obs
# Either pull from row names, index, or from appropriate column in file
obs = my_adata.obs[['cell_id']]

# Store gene-based metadata in a dataframe called var
# Either pull from column names, index, or from appropriate column in file
var = my_adata.var[['gene_id']]
	
# You should now have a matrix called "count_matrix" a sample vector called "obs" and a gene vector called "var"


#######################################################

# PART 2: DATA OUTPUT

# This section describes how to output your variable in a compressed h5ad file, check the size, and then split the output file into multiple files for upload to MapMyCells if the size exceeds 2GB.  All the code below should be run, regardless of which input option was run above.

# Import libraries
import anndata
import os
import scipy.sparse

# Convert count matrix to a CSR (row-based access) sparse matrix and save in anndata format. 		
count_matrix = scipy.sparse.csr_matrix(count_matrix)
ad = anndata.AnnData(
     X=count_matrix,
     obs=obs,
     var=var
)	 

# Write to your location of choice using compression
output_path = 'path/to/file/new_adata_for_MapMyCells.h5ad'

ad.write_h5ad(output_path, compression='gzip')

# Determine and print the file size
file_size_bytes = os.path.getsize(output_path)
print("File size in bytes:", file_size_bytes)
# IF THE PRINTED VALUE IS LESS THAN 2,147,483,648 (2GB), STOP NOW. YOU CAN PROCEED TO MAPMYCELLS UPLOAD!

# If needed, divide your data set into multiple files for upload.  Since the current (November 2024) algorithm are all indenpendently run on each cell, you should get approximately identical answers regardless of how you divide your file.

# Determing the number of files and rows per file
import math
N = math.ceil(file_size_bytes/2147483648)
num_rows = ad.shape[0]
rows_per_subset = num_rows // N
# If N>8 or so, consider using the code-based approach: https://github.com/AllenInstitute/cell_type_mapper

# Write each file to your location of choice using compression
for i in range(N):
     start_idx = i * rows_per_subset
     end_idx = min((i + 1) * rows_per_subset, num_rows)

     ad = anndata.AnnData(
		X=count_matrix[start_idx:end_idx, :],
		obs=obs[start_idx:end_idx],
		var=var
	 )
     ad.write_h5ad(f'path/to/file/new_adata_for_MapMyCells_{i}_.h5ad', compression='gzip')

# You can now upload these files separately to MapMyCells and merge the resulting files together