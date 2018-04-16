import h5py
import numpy as np
import scanpy.api as sc
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_version_and_date()
sc.logging.print_versions_dependencies_numerics()

def remove_cell_cycle(input_file, out_file):

    h5f = h5py.File(input_file, 'r')

    matrix = h5f['matrix'][:]
    gene_names = h5f['gene_attrs']['gene_names'].value
    decoder = np.vectorize(lambda t: t.decode('UTF-8'))
    gene_names = decoder(gene_names)

    adata = sc.AnnData(X=matrix, var=gene_names)

    # Load cell cycle genes defined in [Tirosh et al, 2015](https://doi.org/10.1126/science.aad0501).
    # It is a list of 97 genes, represented by their gene symbol.

    cell_cycle_genes = [x.strip() for x in open('./data/regev_lab_cell_cycle_genes.txt')]

    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var[0].values]

    # this is needed otherwise scanpy cannot tell the index
    adata.var_names = gene_names

    # Log-transformation of data and scaling should always be performed before scoring
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata)

    # calculate the cell cycle scores
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

    sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
    sc.pp.scale(adata)

    matrix = adata.X
    cell_phase = np.array(adata.obs['phase'].values, dtype='S10')

    # write the output
    f = h5py.File(out_file, "w")
    f.create_dataset(name = 'matrix', data = matrix)
    gg = f.create_group('gene_attrs')
    cg = f.create_group('cell_attrs')
    cg.create_dataset(name='cell_phase', data = cell_phase)

    for key in h5f['gene_attrs'].keys():
        gg.create_dataset(name = key, data =h5f['gene_attrs'][key].value)
    for key in h5f['cell_attrs'].keys():
        cg.create_dataset(name = key, data = h5f['cell_attrs'][key].value)

    f.close()
    h5f.close()
