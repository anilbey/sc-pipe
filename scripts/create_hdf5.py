import argparse
import numpy as np
from scipy.io import mmread
import h5py


parser = argparse.ArgumentParser()
parser.add_argument("genes_file", help="tsv file containing gene names")
parser.add_argument("matrix_file", help="file containing the geneXcell matrix")
parser.add_argument("barcodes_file", help="file containing the cell ids(barcodes)")
parser.add_argument("output_file", help="path to the output hdf5 file")
parser.add_argument("--n_threads", help="the number of threads", type=int,
        default=1)
args = parser.parse_args()


def cellranger_to_hdf5 (genes, matrix_file, barcodes, out_path):
    # to read a matrix market file
    matrix = mmread(matrix_file.__str__()).astype("float32").todense().T
    
    gene_names = np.genfromtxt(genes.__str__(), dtype='S16')[:, 1]
    cell_names = np.genfromtxt(barcodes.__str__(), dtype='S16')
   
    f = h5py.File(out_path.__str__(), "w")
    f.create_dataset(name = 'matrix', data = matrix)
    gg = f.create_group('gene_attrs')
    gg.create_dataset(name = 'gene_names', data = gene_names)
    cg = f.create_group('cell_attrs')
    cg.create_dataset(name = 'cell_names', data = cell_names)
    cg.create_dataset(name = 'cells_on_rows', data = True)

    f.close()

cellranger_to_hdf5(args.genes_file, args.matrix_file, args.barcodes_file,
        args.output_file)
