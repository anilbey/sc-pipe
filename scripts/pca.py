from sklearn.decomposition import PCA
import h5py
import numpy as np
import logging
import time

def apply_pca(input_file, output_file, n_components, threads):

    h5f = h5py.File(input_file, 'r')
    matrix = h5f['matrix'][:]
    h5f.close()
    matrix = np.log1p(matrix)
    
    pca = PCA(n_components = n_components)
    pca_Zhat = pca.fit_transform(matrix)
    np.savetxt(output_file, pca_Zhat, delimiter=",")


# try:
apply_pca(snakemake.input.__str__(), snakemake.output.__str__(), snakemake.params.n_components, snakemake.threads)
'''
    except Exception as e:
    log_time = time.strftime("%d.%m.%Y-%H:%M:%S")
    f_name = 'exception/pca'+log_time+'.log'
    logging.basicConfig(filename=f_name,level=logging.INFO)
    logging.exception(e, exc_info=True)
'''