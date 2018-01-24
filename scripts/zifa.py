from ZIFA import ZIFA, block_ZIFA
import h5py
import numpy as np

def apply_zifa(input_file, output_file, n_components, n_blocks, threads):

    h5f = h5py.File(input_file, 'r')
    matrix = h5f['matrix'][:]
    h5f.close()
    matrix = np.log1p(matrix)
   
    print('##ZIFA: number of blocks')
    print(n_blocks)
    Zhat, params = block_ZIFA.fitModel(matrix, n_components, n_blocks = n_blocks)
    np.savetxt(output_file, Zhat, delimiter=",")
    
    
apply_zifa(snakemake.input.__str__(), snakemake.output.__str__(), snakemake.params.n_components, snakemake.params.n_blocks, snakemake.threads)
