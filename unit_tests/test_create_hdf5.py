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
input_path = os.path.join(dirname, '../data/pbmc3k/filtered_gene_bc_matrices/hg19/')
output_path = os.path.join(dirname, '../data/pbmc3k/pbmc3k.hdf5')
logging.debug('out_file_path: ' + output_path)


@pytest.mark.parametrize("genes_file, matrix_file, barcodes_file, out_path", [(input_path + 'genes.tsv', input_path + 'matrix.mtx', input_path + 'barcodes.tsv', output_path)])
@pytest.mark.run(order=1)
def test_cellranger_to_hdf5(genes_file, matrix_file, barcodes_file, out_path):
    cellranger_to_hdf5(genes_file, matrix_file, barcodes_file, out_path)
    # make sure the file is created
    assert os.path.isfile(out_path)
    # make sure the file is not empty
    assert os.stat(out_path).st_size != 0
