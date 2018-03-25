import numpy as np
from scipy.io import mmread
import h5py


def cellranger_to_hdf5 (genes, matrix_file, barcodes, out_path):
    # to read a matrix market file
    matrix = mmread(matrix_file.__str__()).astype("float32").todense().T
    matrix = np.asarray(matrix)
    
    gene_ids = np.genfromtxt(genes.__str__(), dtype='S16')[:,0]
    gene_names = np.genfromtxt(genes.__str__(), dtype='S16')[:, 1]
    cell_names = np.genfromtxt(barcodes.__str__(), dtype='S16')

    # removing all-zero-genes accross all cells
    detected_genes_index = ~(matrix == 0).all(axis=0) 
    matrix = matrix[:,detected_genes_index]
    gene_ids = gene_ids[detected_genes_index]
    gene_names = gene_names[detected_genes_index]


    f = h5py.File(out_path.__str__(), "w")
    f.create_dataset(name = 'matrix', data = matrix)
    gg = f.create_group('gene_attrs')
    gg.create_dataset(name = 'gene_names', data = gene_names)
    gg.create_dataset(name = 'gene_ids', data = gene_ids)
    cg = f.create_group('cell_attrs')
    cg.create_dataset(name = 'cell_names', data = cell_names)
    cg.create_dataset(name = 'cells_on_rows', data = True)

    f.close()
