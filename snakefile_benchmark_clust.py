import h5py
import pandas as pd
import numpy as np
from scipy.io import mmread
import loompy

configfile: 'config/config.json'

SAMPLES = ['melanomaS2']
HDF5_OUTPUT = 'test_hdf5_data'




'''
rules
'''

rule all:
	input:
		simulated_samples = expand(HDF5_OUTPUT+'/{sample}_sim_loc'+'{loc}'+'.loom', sample=SAMPLES, loc=config['splat_simulate']['de_loc_factor'])
	shell:
		'echo test rule all {input.simulated_samples}'
    

rule create_loom:
	input:
		genes_file = expand('cell_ranger_output/{sample}/outs/filtered_gene_bc_matrices/GRCh38/genes.tsv', sample=SAMPLES),
		matrix_file = expand('cell_ranger_output/{sample}/outs/filtered_gene_bc_matrices/GRCh38/matrix.mtx', sample=SAMPLES),
		barcodes_file = expand('cell_ranger_output/{sample}/outs/filtered_gene_bc_matrices/GRCh38/barcodes.tsv', sample=SAMPLES)
	output:
		expand(HDF5_OUTPUT+'/{sample}.loom', sample=SAMPLES)
	script:
		"scripts/create_loom.py"


rule preprocess_zheng17:
	input:
		loom_file = expand(HDF5_OUTPUT+'/{sample}.loom', sample=SAMPLES)	
	output:
		expand(HDF5_OUTPUT+'/{sample}_zheng17.loom', sample=SAMPLES)
	script:
		"scripts/preprocess_zheng17.py"


rule simulate_data:
    input:
        sample_loom = expand(HDF5_OUTPUT+'/{sample}_zheng17.loom', sample=SAMPLES)
    params:
        de_loc_factor = config['splat_simulate']['de_loc_factor']
        #de_prob = config['splat_simulate']['de_prob'],
        #de_dr_prob = config['splat_simulate']['de_dr_prob'],
        #de_scale_factor = config['splat_simulate']['de_scale_factor']	
    output:
        expand(HDF5_OUTPUT+'/{sample}_sim_loc'+'{loc}'+'.loom', sample=SAMPLES, loc=config['splat_simulate']['de_loc_factor'])
    script:
         "scripts/data_simulation.R"
		