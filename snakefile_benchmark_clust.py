

configfile: 'config/config.json'


SAMPLE = 'melanomaS2'

HDF5_OUTPUT = 'test_hdf5_data'
ANALYSIS_OUTPUT = 'analysis'



'''
rules
'''

rule all:
    input:
        expand(ANALYSIS_OUTPUT+'/pca/'+SAMPLE+'_sim_loc'+'{loc}'+'.csv', loc=config['splat_simulate']['de_loc_factor']),
        expand(ANALYSIS_OUTPUT+'/factor_analysis/'+SAMPLE+'_sim_loc'+'{loc}'+'.csv', loc=config['splat_simulate']['de_loc_factor']),
        expand(ANALYSIS_OUTPUT+'/tsne/'+SAMPLE+'_sim_loc'+'{loc}'+'.csv', loc=config['splat_simulate']['de_loc_factor'])
        #expand(ANALYSIS_OUTPUT+'/zifa/'+SAMPLE+'_sim_loc'+'{loc}'+'.csv', loc=config['splat_simulate']['de_loc_factor'])
        
    shell:
        'echo test rule all {input}'

# rule cellranger count (parallel)
# rule cellranger aggr


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
    #params:
        #de_loc_factor = config['splat_simulate']['de_loc_factor']
        #de_prob = config['splat_simulate']['de_prob'],
        #de_dr_prob = config['splat_simulate']['de_dr_prob'],
        #de_scale_factor = config['splat_simulate']['de_scale_factor']	
    output:
       HDF5_OUTPUT+'/'+SAMPLE+'_sim_loc'+'{loc}'+'.loom'
    script:
         "scripts/data_simulation.R"


# run those rules for all simulated data

rule pca:
    input:
        HDF5_OUTPUT+'/'+SAMPLE+'_sim_loc'+'{loc}'+'.loom'
    params:
        n_components = config['dim_reduction']['pca']['n_components']
    output:
        ANALYSIS_OUTPUT+'/pca/'+SAMPLE+'_sim_loc'+'{loc}'+'.csv'
    script:
        "scripts/pca.py"


rule factor_analysis:
    input:
        HDF5_OUTPUT+'/'+SAMPLE+'_sim_loc'+'{loc}'+'.loom'
    params:
        n_components = config['dim_reduction']['factor_analysis']['n_components']
    output:
        ANALYSIS_OUTPUT+'/factor_analysis/'+SAMPLE+'_sim_loc'+'{loc}'+'.csv'
    script:
        "scripts/factor_analysis.py"

rule tsne:
    input:
        HDF5_OUTPUT+'/'+SAMPLE+'_sim_loc'+'{loc}'+'.loom'
    params:
        n_components = config['dim_reduction']['tsne']['n_components']
    output:
        ANALYSIS_OUTPUT+'/tsne/'+SAMPLE+'_sim_loc'+'{loc}'+'.csv'
    script:
        "scripts/tsne.py"

rule zifa:
    input:
        HDF5_OUTPUT+'/'+SAMPLE+'_sim_loc'+'{loc}'+'.loom'
    params:
        n_components = config['dim_reduction']['block_zifa']['n_components'],
        n_blocks = config['dim_reduction']['block_zifa']['n_blocks']
    output:
        ANALYSIS_OUTPUT+'/zifa/'+SAMPLE+'_sim_loc'+'{loc}'+'.csv'
    script:
        "scripts/zifa.py"

'''




rule simlr:


'''