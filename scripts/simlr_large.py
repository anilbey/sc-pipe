import h5py
import time
import numpy as np
import SIMLR
import logging
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="input matrix hdf5 file")
parser.add_argument("output_file", help="path to the output csv file")
parser.add_argument("n_components", help="the number of latent factors for simlr",
        type=int)
parser.add_argument("pca_components", help="# PCs to use", type=int)
parser.add_argument("n_neighbours", help="#nearest neighbours", type=int)
parser.add_argument("max_iter", help="max number of iterations", type=int)
parser.add_argument("--n_threads", help="the number of threads", type=int, default=1)
args = parser.parse_args()

#  https://github.com/anilbey/SIMLR_PY simlr version is used here to support multi-threading
def apply_simlr(input_file, output_file, n_components, pca_components, n_neighbours, max_iter, threads):
    h5f = h5py.File(input_file, 'r')
    matrix = h5f['matrix'][:]
    barcodes = h5f['cell_attrs']['cell_names'].value
    h5f.close()
    
    matrix = np.log1p(matrix)
    X = matrix
    
    # Selecting 500 features with PCA
    if X.shape[1]>500:
    # fast_pca assumes the number of cells > 500 therefore try-catch
        try:
            X = SIMLR.helper.fast_pca(X,pca_components)
        except:
            pass
    
    # Running Simlr 
    simlr = SIMLR.SIMLR_LARGE(num_of_rank=n_components, num_of_neighbor=n_neighbours, max_iter=max_iter)
    S, F,val, ind = simlr.fit(X)
    
    np.savetxt(output_file, F, delimiter=",")
    df = pd.DataFrame(barcodes)
    df = pd.concat([df, pd.DataFrame(F)], axis=1)
    
    df.to_csv(output_file,header=False, index=False)
    

apply_simlr(args.input_file, args.output_file, args.n_components, args.pca_components, args.n_neighbours, args.max_iter, args.n_threads)

    
