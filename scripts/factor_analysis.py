from sklearn.decomposition import FactorAnalysis
import h5py
import numpy as np
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="input matrix hdf5 file")
parser.add_argument("output_file", help="path to the output csv file")
parser.add_argument("n_components", help="the number of latent factors",
        type=int)
parser.add_argument("--n_threads", help="the number of threads", type=int,
        default=1)
args = parser.parse_args()



def apply_fa(input_file, output_file, n_components, threads):

    h5f = h5py.File(input_file, 'r')
    matrix = h5f['matrix'][:]
    barcodes = h5f['cell_attrs']['cell_names'].value
    h5f.close()
    matrix = np.log1p(matrix)
    
    fa = FactorAnalysis(n_components = n_components)
    fa_Zhat = fa.fit_transform(matrix)
    
    df = pd.DataFrame(barcodes)
    df = pd.concat([df, pd.DataFrame(fa_Zhat)], axis=1)
    
    df.to_csv(output_file,header=False, index=False)


apply_fa(args.input_file, args.output_file, args.n_components, args.n_threads)
