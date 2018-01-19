from sklearn.decomposition import FactorAnalysis
import h5py
import numpy as np


def apply_fa(input_file, output_file, n_components, threads):

    h5f = h5py.File(input_file, 'r')
    matrix = h5f['matrix'][:]
    h5f.close()
    matrix = np.log1p(matrix)
    
    fa = FactorAnalysis(n_components = n_components)
    fa_Zhat = fa.fit_transform(matrix)
    np.savetxt(output_file, fa_Zhat, delimiter=",")


apply_fa(snakemake.input.__str__(), snakemake.output.__str__(), snakemake.params.n_components, snakemake.threads)