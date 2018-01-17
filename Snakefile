import h5py
import pandas as pd
import numpy as np
from scipy.io import mmread
import loompy

SAMPLES = ['melanomaS2']
HDF5_OUTPUT = 'test_hdf5_data'




'''
rules
'''

rule all:
	input:
		expanded_samples = expand(HDF5_OUTPUT+'/{sample}.loom', sample=SAMPLES)
	shell:
		'echo test rule all {input.expanded_samples}'
    

rule create_loom:
	input:
		genes_file = expand('cell_ranger_output/{sample}/outs/filtered_gene_bc_matrices/GRCh38/genes.tsv', sample=SAMPLES),
		matrix_file = expand('cell_ranger_output/{sample}/outs/filtered_gene_bc_matrices/GRCh38/matrix.mtx', sample=SAMPLES),
		barcodes_file = expand('cell_ranger_output/{sample}/outs/filtered_gene_bc_matrices/GRCh38/barcodes.tsv', sample=SAMPLES)
	output:
		expand(HDF5_OUTPUT+'/{sample}.loom', sample=SAMPLES)
	script:
		"scripts/create_loom.py"
