import pandas as pd
import numpy as np
from scipy.io import mmread
import loompy
import h5py
from pathlib import Path



def cellranger_to_loom(genes, matrix_file, barcodes, out_path, threads):
    # python code
    #Path(out_path).touch()
    # to read a matrix market file
    matrix = mmread(matrix_file.__str__()).astype("float32").todense()
    gene_names = np.genfromtxt(genes.__str__(), dtype=str)[:, 1]
    cell_names = np.genfromtxt(barcodes.__str__(), dtype=str)
    col_attrs = { "cells": cell_names }
    row_attrs = { "genes": gene_names }
    
    loompy.create(out_path.__str__(), matrix=matrix, row_attrs=row_attrs, col_attrs=col_attrs)
    
cellranger_to_loom(snakemake.input.genes_file, snakemake.input.matrix_file, snakemake.input.barcodes_file, snakemake.output[0], snakemake.threads)