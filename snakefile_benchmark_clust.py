configfile: 'config/config.json'
SAMPLE = ['hgmm100']  #['melanomaS2']
HDF5_OUTPUT = 'hdf5_data'
SIMULATED_DATA_OUTPUT = config['simulated_data_output']
ANALYSIS_OUTPUT = SIMULATED_DATA_OUTPUT+'/analysis'
LOG_FILES = SIMULATED_DATA_OUTPUT+'/log'
CELL_RANGER_OUTPUT_PATH = config['cell_ranger_output']
#e.g. T_CODE=GRCh38, T_CODE=hg19
T_CODE = config['transcriptome_code']

#ruleorder: simulate_data  > preprocess_zheng17
'''
rules
'''

rule all:
    input:
        direct_clustering_results = expand(ANALYSIS_OUTPUT+'/{method}/'+'clusters/{sample}_sim_loc'+'{loc}'+'.csv', loc=config['splat_simulate']['de_loc_factor'], method=config['clustering']['methods_used'], sample=SAMPLE),
        dim_red_results = expand(ANALYSIS_OUTPUT+'/{method}/{sample}_sim_loc{loc}.csv', loc=config['splat_simulate']['de_loc_factor'], method=config['dim_reduction']['methods_used'], sample=SAMPLE),
        dim_red_clustering_results = expand(ANALYSIS_OUTPUT+'/{method}/'+'clusters/{c_method}_{sample}_sim_loc{loc}.csv',
                loc=config['splat_simulate']['de_loc_factor'], method=config['dim_reduction']['methods_used'], c_method=config['dim_reduction']['clustering_methods'], sample=SAMPLE)

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


# silhouette rules run for each dimensionality reduction results
rule silhouette_hierarchical:
    input:
        ANALYSIS_OUTPUT+'/{method}/{sample}_sim_loc'+'{loc}'+'.csv'    
    output:
        ANALYSIS_OUTPUT+'/{method}/clusters/hierarchical_{sample}_sim_loc{loc}.csv'
    log:
        out = LOG_FILES+'/silhouette_hierarchical/{method}_sample_{sample}_loc_{loc}.out',
        err = LOG_FILES+'/silhouette_hierarchical/{method}_sample_{sample}_loc_{loc}.err'
    params:
        affinity = config['clustering']['hierarchical']['affinity'],
        linkage = config['clustering']['hierarchical']['linkage'],
        k_min = config['clustering']['silhouette']['k_min'],
        k_max = config['clustering']['silhouette']['k_max'],
        metric = config['clustering']['silhouette']['metric']
    shell:
        'python scripts/apply_silhouette_hierarchical.py {input} {output} {params.affinity} {params.linkage} {params.k_min} {params.k_max} {params.metric} 2> {log.err} 1> {log.out}'

rule silhouette_kmeans:
    input:
        ANALYSIS_OUTPUT+'/{method}/{sample}_sim_loc'+'{loc}'+'.csv'    
    output:
        ANALYSIS_OUTPUT+'/{method}/clusters/kmeans_{sample}_sim_loc{loc}.csv'
    log:
        out = LOG_FILES+'/silhouette_kmeans/{method}_sample_{sample}_loc_{loc}.out',
        err = LOG_FILES+'/silhouette_kmeans/{method}_sample_{sample}_loc_{loc}.err'
    params:
        n_init = config['clustering']['kmeans']['n_init'],
        n_jobs = config['clustering']['kmeans']['n_jobs'],
        k_min = config['clustering']['silhouette']['k_min'],
        k_max = config['clustering']['silhouette']['k_max'],
        metric = config['clustering']['silhouette']['metric']
    shell:
        'python scripts/apply_silhouette_kmeans.py {input} {output} {params.n_init} {params.n_jobs} {params.k_min} {params.k_max} {params.metric} 2> {log.err} 1> {log.out} '

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
        n_components = config['dim_reduction']['tsne']['n_components'],
        init = config['dim_reduction']['tsne']['init']
    output:
        ANALYSIS_OUTPUT+'/tsne/{sample}_sim_loc'+'{loc}'+'.csv'
    log:
        out = LOG_FILES+'/tsne/sample_{sample}loc_{loc}.out',
        err = LOG_FILES+'/tsne/sample_{sample}loc_{loc}.err'
    shell:
        "python scripts/apply_tsne.py {input} {output} {params.n_components} {params.init}  2> {log.err} 1> {log.out}"

rule phenograph:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_loc{loc}_zheng17.h5'
    params:
        n_neighbours = config['clustering']['phenograph']['n_neighbours']
    output:
        ANALYSIS_OUTPUT+'/phenograph/clusters/{sample}_sim_loc'+'{loc}'+'.csv'
    log:
        out = LOG_FILES+'/phenograph/sample_{sample}loc_{loc}.out',
        err = LOG_FILES+'/phenograph/sample_{sample}loc_{loc}.err'
    threads: 8
    shell:
        "python scripts/apply_phenograph.py {input} {output} {params.n_neighbours} --n_threads {threads} 2> {log.err} 1> {log.out}"

rule griph:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_loc{loc}_zheng17.h5'
    params:
        #TODO read from the config
    output:
        ANALYSIS_OUTPUT+'/griph/clusters/{sample}_sim_loc'+'{loc}'+'.csv'
    log:
        out = LOG_FILES+'/griph/sample_{sample}loc_{loc}.out',
        err = LOG_FILES+'/griph/sample_{sample}loc_{loc}.err'
    threads:8
    shell:
        "Rscript scripts/apply_griph.R --input {input} --output {output} --threads {threads} 2> {log.err} 1> {log.out}"


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





