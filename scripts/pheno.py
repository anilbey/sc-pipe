import phenograph
import h5py
import numpy as np


def apply_pheno(input_file, output_file, n_neighbours, threads):
    h5f = h5py.File(input_file, 'r')
    matrix = h5f['matrix'][:]
    h5f.close()
    matrix = np.log1p(matrix)

    communities, graph, Q = phenograph.cluster(data=matrix,k=n_neighbours,n_jobs=threads)
    np.savetxt(output_file, communities, fmt='%i', delimiter=",")



apply_pheno(snakemake.input.__str__(), snakemake.output.__str__(), snakemake.params.n_neighbours, snakemake.threads)
