from ZIFA import ZIFA, block_ZIFA
import h5py
import numpy as np
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="input matrix hdf5 file")
parser.add_argument("output_file", help="path to the output csv file")
parser.add_argument("n_components", help="the number of principal components",
        type=int)
parser.add_argument("--n_blocks", help="the number of blocks to speedup",
        type=int, default=8)
args = parser.parse_args()


def apply_zifa(input_file, output_file, n_components, n_blocks):

    h5f = h5py.File(input_file, 'r')
    matrix = h5f['matrix'][:]
    barcodes = h5f['cell_attrs']['cell_names'].value
    h5f.close()
    matrix = np.log1p(matrix)
   
    Zhat, params = block_ZIFA.fitModel(matrix, n_components, n_blocks = n_blocks)

    df = pd.DataFrame(barcodes)
    df = pd.concat([df, pd.DataFrame(Zhat)], axis=1)
    df.to_csv(output_file,header=False, index=False)

   
    
apply_zifa(args.input_file, args.output_file, args.n_components, args.n_blocks)
