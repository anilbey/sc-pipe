import h5py
import pandas as pd
import numpy as np
from scipy.io import mmread
import loompy

SAMPLES = ['melanomaS2']
HDF5_OUTPUT = 'test_hdf5_data'


'''
methods
'''
#def cellranger_to_loom(wildcards):
	


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

#	shell:
#		"touch {output}"	
		      
'''
rule create_loom:
	input:
		genes_file = expand('cell_ranger_output/{sample}/outs/filtered_gene_bc_matrices/GRCh38/genes.tsv', sample=SAMPLES),
		matrix_file = expand('cell_ranger_output/{sample}/outs/filtered_gene_bc_matrices/GRCh38/matrix.mtx', sample=SAMPLES),
		barcodes_file = expand('cell_ranger_output/{sample}/outs/filtered_gene_bc_matrices/GRCh38/barcodes.tsv', sample=SAMPLES)
	output:
		expand('test_hdf5_data/{sample}.loom', sample=SAMPLES)
	run:
		# to read a matrix market file
		matrix = mmread(filename_data).astype("float32").todense()
		gene_names = np.genfromtxt(genes, dtype=str)[:, 1]
		cell_names = np.genfromtxt(barcodes, dtype=str)
		# loom row & column conventions
		col_attrs = { "cells": cell_names }
		row_attrs = { "genes": gene_names }
		# output
		loompy.create({output}, matrix=matrix, row_attrs=row_attrs, col_attrs=col_attrs)
'''		