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
from unsupervised.phenograph import Phenograph

# test data I/O
input_file = os.path.join(dirname, '../data/pbmc3k/pbmc3k_preprocessed.hdf5')
output_file = os.path.join(dirname, './test_data/clusters.csv')
n_neighbours = 100
threads = 8

logging.debug('out_file_path: ' + output_file)


@pytest.mark.parametrize("input_file, output_file, n_neighbours, threads", [(input_file, output_file, n_neighbours, threads)])
@pytest.mark.run(order=3)
def test_phenograph(input_file,output_file, n_neighbours, threads):
    pheno = Phenograph(n_neighbours, threads)
    pheno.load_from_hdf5(input_file)
    pheno.log_normalize()
    pheno.apply()
    pheno.write_csv(output_file)
    # make sure the file is created
    assert os.path.isfile(output_file)
    # make sure the file is not empty
    assert os.stat(output_file).st_size != 0
