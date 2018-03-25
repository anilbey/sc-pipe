import pytest
import sys
import os.path
sys.path.append("../scripts")
from utils import cellranger_to_hdf5

# test data paths
pbmc3k = '../data/pbmc3k/filtered_gene_bc_matrices/hg19'
pbmc3k_out = './test_data/pbmc3k.hdf5'


@pytest.mark.parametrize("genes_file, matrix_file, barcodes_file, out_path", [(pbmc3k + 'genes.tsv', pbmc3k + 'matrix.mtx', pbmc3k + 'barcodes.tsv', pbmc3k_out)])
def test_cellranger_to_hdf5(genes_file, matrix_file, barcodes_file, out_path):
    cellranger_to_hdf5(genes_file, matrix_file, barcodes_file, out_path)
    assert os.path.isfile(out_path) 
