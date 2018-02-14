configfile: 'config/config.json'
SAMPLE = ['melanomaS2']
HDF5_OUTPUT = 'hdf5_data'
SIMULATED_DATA_OUTPUT = config['simulated_data_output']
ANALYSIS_OUTPUT = SIMULATED_DATA_OUTPUT+'/analysis'
LOG_FILES = SIMULATED_DATA_OUTPUT+'/log'
CELL_RANGER_OUTPUT_PATH = config['cell_ranger_output']


def serialize(list_var):
#   Serializes the list into a comma separated string
    if isinstance(list_var, list): 
        return (','.join(map(str, list_var)))
    else:
        return list_var


#ruleorder: simulate_data  > preprocess_zheng17
'''
rules
'''

rule all:
    input:
        expand(ANALYSIS_OUTPUT+'/{method}/'+'clusters/{sample}_sim_loc'+'{loc}'+'clusters.csv',
                loc=config['splat_simulate']['de_loc_factor'],method=config['dim_reduction']['methods_used']+config['clustering']['methods_used'],sample=SAMPLE)
    shell:
        'echo test rule all {input}'

# rule cellranger count (parallel)
# rule cellranger aggr

rule cluster_results:
    input:
        ANALYSIS_OUTPUT+'/{method}/{sample}_sim_loc'+'{loc}'+'.csv'
    params:
        n_clusters = len(config['splat_simulate']['group_prob'])
    output:
        ANALYSIS_OUTPUT+'/{method}/'+'clusters/{sample}_sim_loc'+'{loc}'+'clusters.csv'
    log:
        LOG_FILES+'/k-means/sample_{sample}loc_{loc}.log'
    script:
        "scripts/k-means.py"

rule create_hdf5:
    input:
        genes_file = CELL_RANGER_OUTPUT_PATH+'/{sample}/outs/filtered_gene_bc_matrices/GRCh38/genes.tsv',
        matrix_file = CELL_RANGER_OUTPUT_PATH+'/{sample}/outs/filtered_gene_bc_matrices/GRCh38/matrix.mtx',
        barcodes_file = CELL_RANGER_OUTPUT_PATH+'/{sample}/outs/filtered_gene_bc_matrices/GRCh38/barcodes.tsv'
    output:
        HDF5_OUTPUT+'/{sample}.h5'
    log:
        out = LOG_FILES+'/create_hdf5/sample_{sample}.out',
        err = LOG_FILES+'/create_hdf5/sample_{sample}.err'
    shell:
        'python scripts/create_hdf5.py {input.genes_file} {input.matrix_file} {input.barcodes_file} {output} 2> {log.err} 1> {log.out} '

rule simulate_data:
    '''
    group_prob can be a list
    '''
    input:
        sample_hdf5 = rules.create_hdf5.output
    params:
        group_prob = config['splat_simulate']['group_prob'],
        dropout_present = config['splat_simulate']['dropout_present']
    wildcard_constraints:
        loc="\d(\.\d+)?"
    output:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_loc{loc}.h5'
    log:
        LOG_FILES+'/simulate_data/sample_{sample}loc_{loc}.log'
    script:
        "scripts/data_simulation.R"


rule preprocess_zheng17:
    '''
    regex pattern 
        loc: 0.25, 0.5, 1, 1.4, 2
    '''
    input:
        hdf5_file = rules.simulate_data.output 
    params:
        transpose = False
    wildcard_constraints:
        loc="\d(\.\d+)?"
    output:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_loc{loc}_zheng17.h5'
    log:
        LOG_FILES+'/preprocess_zheng17/sample_{sample}loc_{loc}.log'
    script: 
        "scripts/preprocess_zheng17.py"


rule pca:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_loc{loc}_zheng17.h5'
    params:
        n_components = config['dim_reduction']['pca']['n_components']
    output:
        ANALYSIS_OUTPUT+'/pca/{sample}_sim_loc'+'{loc}'+'.csv'
    log:
        out = LOG_FILES+'/pca/sample_{sample}loc_{loc}.out',
        err = LOG_FILES+'/pca/sample_{sample}loc_{loc}.err'
    shell:
        "python scripts/pca.py {input} {output} {params.n_components} 2> {log.err} 1> {log.out}"

rule simlr:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_loc{loc}_zheng17.h5'
    params:
        n_components = config['dim_reduction']['simlr']['n_components'],
        pca_components = config['dim_reduction']['simlr']['pca_components'],
        n_neighbours = config['dim_reduction']['simlr']['n_neighbours'],
        max_iter = config['dim_reduction']['simlr']['max_iter']
    output:
        ANALYSIS_OUTPUT+'/simlr/{sample}_sim_loc'+'{loc}'+'.csv'
    log:
        LOG_FILES+'/simlr/sample_{sample}loc_{loc}.log'
    script:
        "scripts/simlr_large.py"


rule factor_analysis:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_loc{loc}_zheng17.h5'
    params:
        n_components = config['dim_reduction']['factor_analysis']['n_components']
    output:
        ANALYSIS_OUTPUT+'/factor_analysis/{sample}_sim_loc'+'{loc}'+'.csv'
    log:
        LOG_FILES+'/factor_analysis/sample_{sample}loc_{loc}.log' 
    script:
        "scripts/factor_analysis.py"

rule tsne:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_loc{loc}_zheng17.h5'
    params:
        n_components = config['dim_reduction']['tsne']['n_components']
    output:
        ANALYSIS_OUTPUT+'/tsne/{sample}_sim_loc'+'{loc}'+'.csv'
    log:
        LOG_FILES+'/tsne/sample_{sample}loc_{loc}.log'
    script:
        "scripts/tsne.py"

rule phenograph:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_loc{loc}_zheng17.h5'
    params:
        n_neighbours = config['clustering']['phenograph']['n_neighbours']
    output:
        ANALYSIS_OUTPUT+'/phenograph/clusters/{sample}_sim_loc'+'{loc}'+'clusters.csv'
    log:
        LOG_FILES+'/phenograph/sample_{sample}loc_{loc}.log'
    threads:8
    script:
        "scripts/pheno.py"

rule louvain:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_loc{loc}_zheng17.h5'
    params:
        #TODO read from the config
    output:
        ANALYSIS_OUTPUT+'/louvain/clusters/{sample}_sim_loc'+'{loc}'+'clusters.csv'
    log:
        LOG_FILES+'/louvain/sample_{sample}loc_{loc}.log'
    threads:8
    script:
        "scripts/louvain.py"

rule griph:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_loc{loc}_zheng17.h5'
    params:
        #TODO read from the config
    output:
        ANALYSIS_OUTPUT+'/griph/clusters/{sample}_sim_loc'+'{loc}'+'clusters.csv'
    log:
        LOG_FILES+'/griph/sample_{sample}loc_{loc}.log'
    threads:8
    script:
        "scripts/griph.R"




rule zifa:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_loc{loc}_zheng17.h5'
    params:
        n_components = config['dim_reduction']['block_zifa']['n_components'],
        n_blocks = config['dim_reduction']['block_zifa']['n_blocks']
    output:
        ANALYSIS_OUTPUT+'/block_zifa/{sample}_sim_loc'+'{loc}'+'.csv'
    log:
        LOG_FILES+'/zifa/sample_{sample}loc_{loc}.log'
    script:
        "scripts/zifa.py"





