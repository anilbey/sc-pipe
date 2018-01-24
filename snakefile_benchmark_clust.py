configfile: 'config/config.json'
SAMPLE = 'melanomaS2'
HDF5_OUTPUT = 'hdf5_data'
SIMULATED_DATA_OUTPUT = 'simulated/dropout_present'
ANALYSIS_OUTPUT = 'analysis/dropout_present'


'''
rules
'''

rule all:
    input:
        expand(ANALYSIS_OUTPUT+'/{method}/'+'clusters/'+SAMPLE+'_sim_loc'+'{loc}'+'clusters.csv',
                loc=config['splat_simulate']['de_loc_factor'],method=config['dim_reduction']['methods_used'])
    shell:
        'echo test rule all {input}'

# rule cellranger count (parallel)
# rule cellranger aggr

rule cluster_results:
    input:
        ANALYSIS_OUTPUT+'/{method}/'+SAMPLE+'_sim_loc'+'{loc}'+'.csv'
    params:
        n_clusters = len(config['splat_simulate']['group_prob'])
    output:
        ANALYSIS_OUTPUT+'/{method}/'+'clusters/'+SAMPLE+'_sim_loc'+'{loc}'+'clusters.csv'
    script:
        "scripts/k-means.py"

rule create_loom:
    input:
        genes_file = 'cell_ranger_output/'+SAMPLE+'/outs/filtered_gene_bc_matrices/GRCh38/genes.tsv',
        matrix_file = 'cell_ranger_output/'+SAMPLE+'/outs/filtered_gene_bc_matrices/GRCh38/matrix.mtx',
        barcodes_file = 'cell_ranger_output/'+SAMPLE+'/outs/filtered_gene_bc_matrices/GRCh38/barcodes.tsv'
    output:
        HDF5_OUTPUT+'/'+SAMPLE+'.loom'
    script:
        "scripts/create_loom.py"


rule preprocess_zheng17:
    input:
        loom_file = HDF5_OUTPUT+'/'+SAMPLE+'.loom'	
    output:
        HDF5_OUTPUT+'/'+SAMPLE+'_zheng17.loom'
    script: 
        "scripts/preprocess_zheng17.py"

# parallel
rule simulate_data:
    input:
        sample_loom = HDF5_OUTPUT+'/'+SAMPLE+'_zheng17.loom'
    params:
        group_prob = config['splat_simulate']['group_prob'],
        dropout_present = config['splat_simulate']['dropout_present']
    output:
        SIMULATED_DATA_OUTPUT+'/'+SAMPLE+'_sim_loc'+'{loc}'+'.loom'
    script:
        "scripts/data_simulation.R"


# run those rules for all simulated data

rule pca:
    input:
        SIMULATED_DATA_OUTPUT+'/'+SAMPLE+'_sim_loc'+'{loc}'+'.loom'
    params:
        n_components = config['dim_reduction']['pca']['n_components']
    output:
        ANALYSIS_OUTPUT+'/pca/'+SAMPLE+'_sim_loc'+'{loc}'+'.csv'
    script:
        "scripts/pca.py"


rule factor_analysis:
    input:
        SIMULATED_DATA_OUTPUT+'/'+SAMPLE+'_sim_loc'+'{loc}'+'.loom'
    params:
        n_components = config['dim_reduction']['factor_analysis']['n_components']
    output:
        ANALYSIS_OUTPUT+'/factor_analysis/'+SAMPLE+'_sim_loc'+'{loc}'+'.csv'
    script:
        "scripts/factor_analysis.py"

rule tsne:
    input:
        SIMULATED_DATA_OUTPUT+'/'+SAMPLE+'_sim_loc'+'{loc}'+'.loom'
    params:
        n_components = config['dim_reduction']['tsne']['n_components']
    output:
        ANALYSIS_OUTPUT+'/tsne/'+SAMPLE+'_sim_loc'+'{loc}'+'.csv'
    script:
        "scripts/tsne.py"

rule zifa:
    input:
        SIMULATED_DATA_OUTPUT+'/'+SAMPLE+'_sim_loc'+'{loc}'+'.loom'
    params:
        n_components = config['dim_reduction']['block_zifa']['n_components'],
        n_blocks = config['dim_reduction']['block_zifa']['n_blocks']
    output:
        ANALYSIS_OUTPUT+'/block_zifa/'+SAMPLE+'_sim_loc'+'{loc}'+'.csv'
    script:
        "scripts/zifa.py"

'''




rule simlr:


'''
