PBMC_PATH = 'data/pbmc68k/filtered_matrices_mex/hg19'
HDF5_OUTPUT = 'hdf5_data'
configfile: 'config/config.json'
SAMPLE = 'PBMC68k'
ANALYSIS_OUTPUT = 'analysis/'+SAMPLE


rule all:
    input:
        expand(ANALYSIS_OUTPUT+'/{method}_clusters.csv',
                method=config['dim_reduction']['methods_used'])
    shell:
        'echo input for rule all: {input}'

rule create_loom:
    input:
        genes_file = PBMC_PATH +'/genes.tsv',
        matrix_file = PBMC_PATH +'/matrix.mtx',
        barcodes_file = PBMC_PATH +'/barcodes.tsv'
    output:
        HDF5_OUTPUT+'/'+SAMPLE+'.loom'
    script:
        "scripts/create_loom.py"

rule cluster_results:
    input:
        ANALYSIS_OUTPUT+'/{method}.csv'
    params:
        n_clusters = 11 # len(config['splat_simulate']['group_prob'])
    output:
        ANALYSIS_OUTPUT+'/{method}_clusters.csv'
    script:
        "scripts/k-means.py"


rule preprocess_zheng17:
    input:
        loom_file = HDF5_OUTPUT+'/'+SAMPLE+'.loom'	
    params:
        transpose = True
    output:
        HDF5_OUTPUT+'/'+SAMPLE+'_zheng17.loom'
    script: 
        "scripts/preprocess_zheng17.py"

rule pca:
    input:
        HDF5_OUTPUT+'/'+SAMPLE+'_zheng17.loom'
    params:
        n_components = config['dim_reduction']['pca']['n_components']
    output:
        ANALYSIS_OUTPUT+'/pca.csv'
    script:
        "scripts/pca.py"

rule simlr:
    input:
        HDF5_OUTPUT+'/'+SAMPLE+'_zheng17.loom'
    params:
        n_components = config['dim_reduction']['simlr']['n_components'],
        pca_components = config['dim_reduction']['simlr']['pca_components'],
        n_neighbours = config['dim_reduction']['simlr']['n_neighbours'],
        max_iter = config['dim_reduction']['simlr']['max_iter']
    output:
        ANALYSIS_OUTPUT+'/simlr.csv'
    script:
        "scripts/simlr_large.py"

rule phenograph:
    input:
        HDF5_OUTPUT+'/'+SAMPLE+'_zheng17.loom'
    params:
        n_neighbours = config['clustering']['phenograph']['n_neighbours']
    output:
        ANALYSIS_OUTPUT+'/phenograph_clusters.csv'
    threads:8
    script:
        "scripts/pheno.py"

rule factor_analysis:
    input:
        HDF5_OUTPUT+'/'+SAMPLE+'_zheng17.loom'
    params:
        n_components = config['dim_reduction']['factor_analysis']['n_components']
    output:
        ANALYSIS_OUTPUT+'/factor_analysis.csv'
    script:
        "scripts/factor_analysis.py"

rule tsne:
    input:
        HDF5_OUTPUT+'/'+SAMPLE+'_zheng17.loom'
    params:
        n_components = config['dim_reduction']['tsne']['n_components']
    output:
        ANALYSIS_OUTPUT+'/tsne.csv'
    script:
        "scripts/tsne.py"

rule zifa:
    input:
        HDF5_OUTPUT+'/'+SAMPLE+'_zheng17.loom'
    params:
        n_components = config['dim_reduction']['block_zifa']['n_components'],
        n_blocks = config['dim_reduction']['block_zifa']['n_blocks']
    output:
        ANALYSIS_OUTPUT+'/block_zifa.csv'
    script:
        "scripts/zifa.py"



