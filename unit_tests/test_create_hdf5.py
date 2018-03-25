import pytest
import sys
import os
import logging
logging.basicConfig(level=logging.DEBUG)

dirname = os.path.dirname(__file__)
logging.debug('__file__ variable: ' +  __file__)
logging.debug('dirname: ' + dirname)
scripts_path = os.path.join(dirname, '../scripts')
sys.path.append(scripts_path)
from utils import cellranger_to_hdf5

# test data paths
pbmc3k = os.path.join(dirname, '../data/pbmc3k/filtered_gene_bc_matrices/hg19')
pbmc3k_out = os.path.join(dirname, './test_data/pbmc3k.hdf5')


@pytest.mark.parametrize("genes_file, matrix_file, barcodes_file, out_path", [(pbmc3k + 'genes.tsv', pbmc3k + 'matrix.mtx', pbmc3k + 'barcodes.tsv', pbmc3k_out)])
def test_cellranger_to_hdf5(genes_file, matrix_file, barcodes_file, out_path):
    cellranger_to_hdf5(genes_file, matrix_file, barcodes_file, out_path)
    assert os.path.isfile(out_path)

