import pandas as pd
import numpy as np
from scipy.io import mmread
import h5py



def cellranger_to_hdf5(genes, matrix_file, barcodes, out_path, threads):
    # to read a matrix market file
    matrix = mmread(matrix_file.__str__()).astype("float32").todense().T
    gene_names = np.genfromtxt(genes.__str__(), dtype='S16')[:, 1]
    cell_names = np.genfromtxt(barcodes.__str__(), dtype='S16')
   
    f = h5py.File(out_path.__str__(), "w")
    f.create_dataset(name = 'matrix', data = matrix)
    f.create_dataset(name = 'gene_names', data = gene_names)
    f.create_dataset(name = 'cell_names', data = cell_names)
    f.create_dataset(name = 'cells_on_rows', data = True)

    f.close()
    
cellranger_to_hdf5(snakemake.input.genes_file, snakemake.input.matrix_file, snakemake.input.barcodes_file, snakemake.output, snakemake.threads)
