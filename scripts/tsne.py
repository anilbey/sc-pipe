from sklearn.manifold import TSNE
import h5py
import numpy as np


def apply_tsne(input_file, output_file, n_components, threads):

    h5f = h5py.File(input_file, 'r')
    matrix = h5f['matrix'][:]
    h5f.close()
    matrix = np.log1p(matrix)
    
    tsne_Zhat = TSNE(init='pca', n_components=n_components).fit_transform(matrix)
    np.savetxt(output_file, tsne_Zhat, delimiter=",")


apply_tsne(snakemake.input.__str__(), snakemake.output.__str__(), snakemake.params.n_components, snakemake.threads)