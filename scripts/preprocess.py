import h5py
import scanpy.api as sc

def preprocess(hdf5_file, out_path, n_top_genes):

    h5f = h5py.File(hdf5_file,'r')
    matrix = h5f['matrix'].value

    adata = sc.AnnData(matrix)

    print(adata.X.shape)
    # sc.pp.normalize_per_cell(adata)          # normalize with total UMI count per cell
    print(adata.X.shape)
    filter_result = sc.pp.filter_genes_dispersion(adata.X, flavor='cell_ranger', n_top_genes=n_top_genes, log=False)
    # filter results is a recarray
    # mask2 to select the top 1000 genes
    mask2 = filter_result.gene_subset
    adata = adata[:, mask2]
    # sc.pp.normalize_per_cell(adata)  # need to redo normalization after filtering

    # Writing the output hdf5 files
    matrix = adata.X

    f = h5py.File(out_path, "w")
    f.create_dataset(name = 'matrix', data = matrix)
    gg = f.create_group('gene_attrs')
    cg = f.create_group('cell_attrs')
    print(h5f['gene_attrs'].keys())
    for key in h5f['gene_attrs'].keys():
        # apply the masks to the gene attributes
        gg.create_dataset(name = key, data =h5f['gene_attrs'][key].value[mask2])
    for key in h5f['cell_attrs'].keys():
        cg.create_dataset(name = key, data = h5f['cell_attrs'][key].value)

    f.close()
    h5f.close()
