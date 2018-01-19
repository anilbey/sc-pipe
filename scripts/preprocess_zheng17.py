import loompy
import scanpy.api as sc



def preprocess_zheng17(loom_file, out_path, threads):
    ds = loompy.connect(loom_file)
    # Scanpy filtering functions assume that: genes=variables, Transpose the matrix
    adata = sc.AnnData(ds[:,:].T, obs=ds.col_attrs, var=ds.row_attrs)
    ds.close()
    sc.pp.filter_genes(adata, min_counts=1)  # only consider genes with more than 1 count
    sc.pp.normalize_per_cell(adata)          # normalize with total UMI count per cell
    filter_result = sc.pp.filter_genes_dispersion(adata.X, flavor='cell_ranger', n_top_genes=1000, log=False)
    adata = adata[:, filter_result.gene_subset]

    # Writing the output loom file
    matrix = adata.X.T
    gene_names = adata.var['genes'].values
    cell_names = adata.obs['cells'].values
    col_attrs = { "cells": cell_names }
    row_attrs = { "genes": gene_names }
    loompy.create(out_path, matrix=matrix, row_attrs=row_attrs, col_attrs=col_attrs)
    
preprocess_zheng17(snakemake.input.loom_file.__str__(), snakemake.output.__str__(), snakemake.threads)