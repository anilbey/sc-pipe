import h5py
import scanpy.api as sc



def preprocess_zheng17(hdf5_file, out_path, threads):
    
    h5f = h5py.File(hdf5_file,'r')

    # Scanpy filtering functions assume that: genes=variables, Transpose the matrix
    adata = sc.AnnData(h5f['matrix'].value, obs=h5f['cell_groups'].value,
            var=h5f['gene_names'].value)
    h5f.close()
    sc.pp.filter_genes(adata, min_counts=1)  # only consider genes with more than 1 count
    sc.pp.normalize_per_cell(adata)          # normalize with total UMI count per cell
    filter_result = sc.pp.filter_genes_dispersion(adata.X, flavor='cell_ranger', n_top_genes=1000, log=False)
    adata = adata[:, filter_result.gene_subset]
    sc.pp.normalize_per_cell(adata)  # need to redo normalization after filtering 

    # Writing the output hdf5 file
    gene_names = adata.var
    cell_groups = adata.obs
    matrix = adata.X
    
    f = h5py.File(out_path, "w")
    f.create_dataset(name = 'matrix', data = matrix)
    # use gene_names[0] because scanpy creates a dataframe with additional
    # informative columns that are not needed 
    f.create_dataset(name = 'gene_names', data = gene_names[0].values)
    f.create_dataset(name = 'cell_groups', data = cell_groups[0].values)

    f.close()
 
    
preprocess_zheng17(snakemake.input.hdf5_file.__str__(),
        snakemake.output.__str__(), snakemake.threads)
