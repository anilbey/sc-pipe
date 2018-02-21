configfile: 'config/config.json'
SAMPLE = ['hgmm100']  #['melanomaS2']
HDF5_OUTPUT = 'hdf5_data'
SIMULATED_DATA_OUTPUT = config['simulated_data_output']
ANALYSIS_OUTPUT = SIMULATED_DATA_OUTPUT+'/analysis'
LOG_FILES = SIMULATED_DATA_OUTPUT+'/log'
CELL_RANGER_OUTPUT_PATH = config['cell_ranger_output']
#e.g. T_CODE=GRCh38, T_CODE=hg19
T_CODE = config['transcriptome_code']

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
        #config['cell_ranger_output']
        expand(ANALYSIS_OUTPUT+'/{method}/'+'clusters/{sample}_sim_loc'+'{loc}'+'clusters.csv', loc=config['splat_simulate']['de_loc_factor'],method=config['dim_reduction']['methods_used']+config['clustering']['methods_used'],sample=SAMPLE)
    shell:
        'echo test rule all {input}'

# run it for each library them perform aggr
rule cellranger_count: # (parallel)
    input:
        fastqs_dir = config['input_fastqs'],
        reference = config['reference_transcriptome'],
        id = config['unique_run_id']
    output:
        config['cell_ranger_output']
    log:
        out = LOG_FILES+'/cellranger_count/sample_{sample}.out',
        err = LOG_FILES+'/cellranger_count/sample_{sample}.err'
    shell:
        'cellranger count --id={input.id} --transcriptome={input.reference} --fastqs={input.fastqs_dir} --nosecondary 2> {log.err} 1> {log.out}'

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
        genes_file =
        CELL_RANGER_OUTPUT_PATH+'/{sample}/outs/filtered_gene_bc_matrices/'+T_CODE+'/genes.tsv',
        matrix_file = CELL_RANGER_OUTPUT_PATH+'/{sample}/outs/filtered_gene_bc_matrices/'+T_CODE+'/matrix.mtx',
        barcodes_file = CELL_RANGER_OUTPUT_PATH+'/{sample}/outs/filtered_gene_bc_matrices/'+T_CODE+'/barcodes.tsv'
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
        out = LOG_FILES+'/simulate_data/sample_{sample}loc_{loc}.out',
        err = LOG_FILES+'/simulate_data/sample_{sample}loc_{loc}.err'
    shell:
        "Rscript scripts/data_simulation.R --input {input.sample_hdf5} --group_prob {params.group_prob} --dropout_present {params.dropout_present} --output {output} --loc {wildcards.loc} 2> {log.err} 1> {log.out}"


rule preprocess_zheng17:
    '''
    regex pattern 
        loc: 0.25, 0.5, 1, 1.4, 2
    '''
    input:
        hdf5_file = rules.simulate_data.output 
    wildcard_constraints:
        loc="\d(\.\d+)?"
    output:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_loc{loc}_zheng17.h5'
    log:
        out = LOG_FILES+'/preprocess_zheng17/sample_{sample}loc_{loc}.out',
        err = LOG_FILES+'/preprocess_zheng17/sample_{sample}loc_{loc}.err'
    shell: 
        "python scripts/preprocess_zheng17.py {input.hdf5_file} {output} 2> {log.err} 1> {log.out}"


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
        "python scripts/apply_pca.py {input} {output} {params.n_components} 2> {log.err} 1> {log.out}"

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
        out = LOG_FILES+'/simlr/sample_{sample}loc_{loc}.out',
        err = LOG_FILES+'/simlr/sample_{sample}loc_{loc}.err'
    shell:
        "python scripts/apply_simlr.py {input} {output} {params.n_components} {params.pca_components} {params.n_neighbours} {params.max_iter} 2> {log.err} 1> {log.out}"


rule factor_analysis:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_loc{loc}_zheng17.h5'
    params:
        n_components = config['dim_reduction']['factor_analysis']['n_components']
    output:
        ANALYSIS_OUTPUT+'/factor_analysis/{sample}_sim_loc'+'{loc}'+'.csv'
    log:
        out = LOG_FILES+'/factor_analysis/sample_{sample}loc_{loc}.out',
        err = LOG_FILES+'/factor_analysis/sample_{sample}loc_{loc}.err'
    shell:
        "python scripts/apply_fa.py {input} {output} {params.n_components} 2> {log.err} 1> {log.out} "

rule tsne:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_loc{loc}_zheng17.h5'
    params:
        n_components = config['dim_reduction']['tsne']['n_components']
    output:
        ANALYSIS_OUTPUT+'/tsne/{sample}_sim_loc'+'{loc}'+'.csv'
    log:
        out = LOG_FILES+'/tsne/sample_{sample}loc_{loc}.out',
        err = LOG_FILES+'/tsne/sample_{sample}loc_{loc}.err'
    shell:
        "python scripts/apply_tsne.py {input} {output} {params.n_components} 2> {log.err} 1> {log.out}"

rule phenograph:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_loc{loc}_zheng17.h5'
    params:
        n_neighbours = config['clustering']['phenograph']['n_neighbours']
    output:
        ANALYSIS_OUTPUT+'/phenograph/clusters/{sample}_sim_loc'+'{loc}'+'clusters.csv'
    log:
        out = LOG_FILES+'/phenograph/sample_{sample}loc_{loc}.out',
        err = LOG_FILES+'/phenograph/sample_{sample}loc_{loc}.err'
    threads: 8
    shell:
        "python scripts/apply_phenograph.py {input} {output} {params.n_neighbours} --n_threads {threads} 2> {log.err} 1> {log.out}"
'''
rule griph:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_loc{loc}_zheng17.h5'
    params:
        #TODO read from the config
    output:
        ANALYSIS_OUTPUT+'/griph/clusters/{sample}_sim_loc'+'{loc}'+'clusters.csv'
    log:
        out = LOG_FILES+'/griph/sample_{sample}loc_{loc}.out',
        err = LOG_FILES+'/griph/sample_{sample}loc_{loc}.err'
    threads:8
    shell:
        "Rscript scripts/griph.R"
'''

rule zifa:
    # by default n_blocks are estimated as n_cells/500 
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_loc{loc}_zheng17.h5'
    params:
        n_components = config['dim_reduction']['block_zifa']['n_components'],
        n_blocks = config['dim_reduction']['block_zifa']['n_blocks']
    output:
        ANALYSIS_OUTPUT+'/block_zifa/{sample}_sim_loc'+'{loc}'+'.csv'
    log:
        out = LOG_FILES+'/zifa/sample_{sample}loc_{loc}.out',
        err = LOG_FILES+'/zifa/sample_{sample}loc_{loc}.err'
    shell:
        "python scripts/apply_zifa.py {input} {output} {params.n_components} --n_blocks {params.n_blocks} 1> {log.out} 2> {log.err}"





