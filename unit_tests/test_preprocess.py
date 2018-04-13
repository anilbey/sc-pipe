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
from preprocess import preprocess

# test data I/O
input_file = os.path.join(dirname, '../data/pbmc3k/pbmc3k.hdf5')
output_file = os.path.join(dirname, '../data/pbmc3k/pbmc3k_preprocessed.hdf5')
n_top_genes = 1000
logging.debug('out_file_path: ' + output_file)


@pytest.mark.parametrize("hdf5_file, out_path, n_top_genes", [(input_file, output_file, n_top_genes)])
@pytest.mark.run(order=2)
def test_preprocess(hdf5_file, out_path, n_top_genes):
    preprocess(hdf5_file, out_path, n_top_genes)
    # make sure the file is created
    assert os.path.isfile(out_path)
    # make sure the file is not empty
    assert os.stat(out_path).st_size != 0
