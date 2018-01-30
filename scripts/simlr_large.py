import h5py
import time
import numpy as np
import SIMLR
import logging
import os

def apply_simlr(input_file, output_file, n_components, pca_components, n_neighbours, max_iter, threads):
    h5f = h5py.File(input_file, 'r')
    matrix = h5f['matrix'][:]
    h5f.close()
    
    matrix = np.log1p(matrix)
    X = matrix
    
    # Selecting 500 features with PCA
    if X.shape[1]>500:
        X = SIMLR.helper.fast_pca(X,pca_components)
    else:
        X = X.todense()
    
    # Running Simlr 
    simlr = SIMLR.SIMLR_LARGE(num_of_rank=n_components, num_of_neighbor=n_neighbours, max_iter=max_iter)
    S, F,val, ind = simlr.fit(X)
    
    np.savetxt(output_file, F, delimiter=",")
    

try:
    apply_simlr(snakemake.input.__str__(), snakemake.output.__str__(), snakemake.params.n_components, snakemake.params.pca_components, snakemake.params.n_neighbours, snakemake.params.max_iter, snakemake.threads)

except Exception as e:
    log_time = time.strftime("%d.%m.%Y-%H:%M:%S")
    f_name = 'exception/simlr'+log_time+'.log'
    logging.basicConfig(filename=f_name,level=logging.INFO)
    logging.exception(e, exc_info=True)

    
