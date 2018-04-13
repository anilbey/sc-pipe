import glob
import numpy as np

configfile: 'config/benchmarking_pipeline_config.json'
SAMPLE = ['melanomaS2']  #['hgmm100']
HDF5_OUTPUT = 'intermediate_files'
SIMULATED_DATA_OUTPUT = config['simulated_data_output']
ANALYSIS_OUTPUT = SIMULATED_DATA_OUTPUT+'/analysis'
LOG_FILES = SIMULATED_DATA_OUTPUT+'/log'
CELL_RANGER_OUTPUT_PATH = config['cell_ranger_output']
#e.g. T_CODE=GRCh38, T_CODE=hg19
T_CODE = config['transcriptome_code']
RUN_ID = config['unique_run_id']
#ruleorder: simulate_data  > preprocess

dim_red_methods_with_fixed_params = config['dim_reduction']['methods_used']
dim_red_methods_with_fixed_params.remove('pca')

clustering_methods_with_fixed_params = config['clustering']['methods_used']
clustering_methods_with_fixed_params.remove('phenograph')

'''
rules
'''

def dirichlet_group_prob(size):
    res = None
    while True:
        dirr = np.random.dirichlet(np.ones(size)*10, size=1)
        if np.sum(dirr)==1.0:
            res = dirr.flatten()
            break
        else:
            continue
    assert np.sum(res)==1.0
    return res.tolist()


def get_all_fastqs(path):
    fastqs =  glob.glob(path+'/*.fastq.gz')

rule all:
    input:
        direct_clustering_results = expand(ANALYSIS_OUTPUT+'/{method}/'+'clusters/{sample}_sim_de{de_prob}_loc'+'{loc}'+'.csv',de_prob= config['splat_simulate']['de_prob'], \
                loc=config['splat_simulate']['de_loc_factor'], method=clustering_methods_with_fixed_params, sample=SAMPLE),

        dim_red_results = expand(ANALYSIS_OUTPUT+'/{method}/{sample}_sim_de{de_prob}_loc{loc}.csv', de_prob=config['splat_simulate']['de_prob'], \
                loc=config['splat_simulate']['de_loc_factor'], method=dim_red_methods_with_fixed_params, sample=SAMPLE),

        phenograph_variations = expand(ANALYSIS_OUTPUT+'/phenograph_{nn}/'+'clusters/{sample}_sim_de{de_prob}_loc'+'{loc}'+'.csv',de_prob= config['splat_simulate']['de_prob'], \
                loc=config['splat_simulate']['de_loc_factor'], nn = config['clustering']['phenograph']['benchmark_n_neighbours'], sample=SAMPLE),

        pca_variations =expand(ANALYSIS_OUTPUT+'/silhouette-pca_{pca_comps}/'+'clusters/{c_method}_{sample}_sim_de{de_prob}_loc{loc}.csv',de_prob=config['splat_simulate']['de_prob'], \
                loc=config['splat_simulate']['de_loc_factor'], pca_comps=config['dim_reduction']['pca']['benchmark_n_components'],c_method=config['dim_reduction']['clustering_methods'],sample=SAMPLE),

        dim_red_clustering_results = expand(ANALYSIS_OUTPUT+'/silhouette-{method}/'+'clusters/{c_method}_{sample}_sim_de{de_prob}_loc{loc}.csv', de_prob=config['splat_simulate']['de_prob'], \
                loc=config['splat_simulate']['de_loc_factor'], method=dim_red_methods_with_fixed_params, c_method=config['dim_reduction']['clustering_methods'], sample=SAMPLE)


#TODO input function        fastqs = get_all_fastqs(config['input_fastqs'])
    shell:
        'echo test rule all {input}'


# run it for each library them perform aggr
# TODO: use wildcard {sample}, unique_run_id must come from {sample}
# TODO: bind it with the create_hdf5 rule, create an output dependency

'''
rule cellranger_count: # (parallel)
    input:
        fastqs_dir = config['input_fastqs']+'/{sample}',
        reference = config['reference_transcriptome']
    params:
        id = config['unique_run_id']
    output:
        genes_file = CELL_RANGER_OUTPUT_PATH+'/{sample}/outs/filtered_gene_bc_matrices/'+T_CODE+'/genes.tsv',
        matrix_file = CELL_RANGER_OUTPUT_PATH+'/{sample}/outs/filtered_gene_bc_matrices/'+T_CODE+'/matrix.mtx',
        barcodes_file = CELL_RANGER_OUTPUT_PATH+'/{sample}/outs/filtered_gene_bc_matrices/'+T_CODE+'/barcodes.tsv'
    log:
        out = LOG_FILES+'/cellranger_count/sample_{sample}.out',
        err = LOG_FILES+'/cellranger_count/sample_{sample}.err'
    shell:
        'cellranger count --id={input.id} --transcriptome={input.reference} --fastqs={input.fastqs_dir} --nosecondary 2> {log.err} 1> {log.out}'
'''

# rule cellranger aggr


rule create_hdf5:
    input:
        #genes_file = rules.cellranger_count.output.genes_file,
        genes_file = CELL_RANGER_OUTPUT_PATH+'/'+RUN_ID+'/outs/filtered_gene_bc_matrices/'+T_CODE+'/genes.tsv',
        #matrix_file = rules.cellranger_count.output.matrix_file,
        matrix_file = CELL_RANGER_OUTPUT_PATH+'/'+RUN_ID+'/outs/filtered_gene_bc_matrices/'+T_CODE+'/matrix.mtx',
        #barcodes_file = rules.cellranger_count.output.barcodes_file
        barcodes_file = CELL_RANGER_OUTPUT_PATH+'/'+RUN_ID+'/outs/filtered_gene_bc_matrices/'+T_CODE+'/barcodes.tsv'
    output:
        HDF5_OUTPUT+'/raw_{sample}.h5'
    log:
        out = LOG_FILES+'/create_hdf5/sample_{sample}.out',
        err = LOG_FILES+'/create_hdf5/sample_{sample}.err'
    shell:
        'python scripts/create_hdf5.py -g {input.genes_file} -m {input.matrix_file} -b {input.barcodes_file} -o {output} 2> {log.err} 1> {log.out} '

rule estimate_params:
    input:
        sample_hdf5 = HDF5_OUTPUT+'/raw_{sample}.h5'
    output:
        HDF5_OUTPUT+'/{sample}_estimated.rda'
    log:
        out = LOG_FILES+'/estimate_params/sample_{sample}.out',
        err = LOG_FILES+'/estimate_params/sample_{sample}.err'
    shell:
        'Rscript scripts/estimate_params.R --input {input.sample_hdf5} --output {output} 2> {log.err} 1> {log.out}'

rule simulate_data:
    '''
    group_prob can be a list
    '''
    input:
        sample_hdf5 = ancient(rules.estimate_params.output)
    params:
        group_prob = dirichlet_group_prob(config['splat_simulate']['group_sizes']),
        dropout_present = config['splat_simulate']['dropout_present']
    wildcard_constraints:
        loc="\d(\.\d+)?",
        de_prob="\d(\.\d+)?"
    output:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_de{de_prob}_loc{loc}.h5'
    log:
        out = LOG_FILES+'/simulate_data/sample_{sample}_de{de_prob}_loc_{loc}.out',
        err = LOG_FILES+'/simulate_data/sample_{sample}_de{de_prob}_loc_{loc}.err'
    shell:
        "Rscript scripts/data_simulation.R --input {input.sample_hdf5} --group_prob {params.group_prob} --dropout_present \
         {params.dropout_present} --output {output} --loc {wildcards.loc} --de_prob {wildcards.de_prob} 2> {log.err} 1> {log.out}"


rule preprocess:
    '''
    regex pattern
        loc: 0.25, 0.5, 1, 1.4, 2
    '''
    input:
        hdf5_file = rules.simulate_data.output
    wildcard_constraints:
        loc="\d(\.\d+)?"
    params:
        n_top_genes = config['preprocess']['zheng17']['n_top_genes']
    output:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_de{de_prob}_loc{loc}_zheng17.h5'
    log:
        out = LOG_FILES+'/preprocess/sample_{sample}_de{de_prob}_loc_{loc}.out',
        err = LOG_FILES+'/preprocess/sample_{sample}_de{de_prob}_loc_{loc}.err'
    shell:
        "python scripts/apply_preprocess.py -i {input.hdf5_file} -o {output} --n_top_genes {params.n_top_genes} 2> {log.err} 1> {log.out}"

# silhouette rules run for each dimensionality reduction results
rule silhouette_hierarchical:
    input:
        ANALYSIS_OUTPUT+'/{method}/{sample}_sim_de{de_prob}_loc'+'{loc}'+'.csv'
    output:
        ANALYSIS_OUTPUT+'/silhouette-{method}/clusters/hierarchical_{sample}_sim_de{de_prob}_loc{loc}.csv'
    log:
        out = LOG_FILES+'/silhouette_hierarchical/{method}_sample_{sample}_de{de_prob}_loc_{loc}.out',
        err = LOG_FILES+'/silhouette_hierarchical/{method}_sample_{sample}_de{de_prob}_loc_{loc}.err'
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
        ANALYSIS_OUTPUT+'/{method}/{sample}_sim_de{de_prob}_loc'+'{loc}'+'.csv'
    output:
        ANALYSIS_OUTPUT+'/silhouette-{method}/clusters/kmeans_{sample}_sim_de{de_prob}_loc{loc}.csv'
    log:
        out = LOG_FILES+'/silhouette_kmeans/{method}_sample_{sample}_de{de_prob}_loc_{loc}.out',
        err = LOG_FILES+'/silhouette_kmeans/{method}_sample_{sample}_de{de_prob}_loc_{loc}.err'
    params:
        n_init = config['clustering']['kmeans']['n_init'],
        n_jobs = config['clustering']['kmeans']['n_jobs'],
        k_min = config['clustering']['silhouette']['k_min'],
        k_max = config['clustering']['silhouette']['k_max'],
        metric = config['clustering']['silhouette']['metric']
    shell:
        'python scripts/apply_silhouette_kmeans.py {input} {output} {params.n_init} {params.n_jobs} {params.k_min} {params.k_max} {params.metric} 2> {log.err} 1> {log.out} '

rule pca:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_de{de_prob}_loc{loc}_zheng17.h5'
    params:
        # this information is already contained in the wildcard
        # n_components = config['dim_reduction']['pca']['n_components']
    output:
        ANALYSIS_OUTPUT+'/pca_{pca_comps}/{sample}_sim_de{de_prob}_loc'+'{loc}'+'.csv'
    log:
        out = LOG_FILES+'/pca/sample_{sample}_de{de_prob}loc_{loc}.out',
        err = LOG_FILES+'/pca/sample_{sample}_de{de_prob}loc_{loc}.err'
    shell:
        "python scripts/apply_pca.py {input} {output} {wildcards.pca_comps} 2> {log.err} 1> {log.out}"

rule phenograph:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_de{de_prob}_loc{loc}_zheng17.h5'
    params:
        #n_neighbours = config['clustering']['phenograph']['n_neighbours']
    output:
        ANALYSIS_OUTPUT+'/phenograph_{nn}/clusters/{sample}_sim_de{de_prob}_loc'+'{loc}'+'.csv'
    log:
        out = LOG_FILES+'/phenograph/sample_{sample}_de{de_prob}loc_{loc}.out',
        err = LOG_FILES+'/phenograph/sample_{sample}_de{de_prob}loc_{loc}.err'
    threads: 8
    shell:
        "python scripts/apply_phenograph.py {input} {output} {wildcards.nn} --n_threads {threads} 2> {log.err} 1> {log.out}"

rule simlr:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_de{de_prob}_loc{loc}_zheng17.h5'
    params:
        n_components = config['clustering']['simlr']['n_components'],
        pca_components = config['clustering']['simlr']['pca_components'],
        n_neighbours = config['clustering']['simlr']['n_neighbours'],
        max_iter = config['clustering']['simlr']['max_iter'],
        # silhouette parameters
        k_min = config['clustering']['silhouette']['k_min'],
        k_max = config['clustering']['silhouette']['k_max'],
        metric = config['clustering']['silhouette']['metric']
    output:
        ANALYSIS_OUTPUT+'/simlr/clusters/{sample}_sim_de{de_prob}_loc'+'{loc}'+'.csv'
    log:
        out = LOG_FILES+'/simlr/sample_{sample}_de{de_prob}loc_{loc}.out',
        err = LOG_FILES+'/simlr/sample_{sample}_de{de_prob}loc_{loc}.err'
    shell:
        "python scripts/apply_silhouette_simlr.py --input {input} --output {output} --n_components {params.n_components} --pca_components {params.pca_components} \
         --n_neighbours {params.n_neighbours} --max_iter {params.max_iter} --k_min {params.k_min} --k_max {params.k_max} --metric {params.metric}  2> {log.err} 1> {log.out}"


rule factor_analysis:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_de{de_prob}_loc{loc}_zheng17.h5'
    params:
        n_components = config['dim_reduction']['factor_analysis']['n_components']
    output:
        ANALYSIS_OUTPUT+'/factor_analysis/{sample}_sim_de{de_prob}_loc'+'{loc}'+'.csv'
    log:
        out = LOG_FILES+'/factor_analysis/sample_{sample}_de{de_prob}loc_{loc}.out',
        err = LOG_FILES+'/factor_analysis/sample_{sample}_de{de_prob}loc_{loc}.err'
    shell:
        "python scripts/apply_fa.py {input} {output} {params.n_components} 2> {log.err} 1> {log.out} "

rule tsne:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_de{de_prob}_loc{loc}_zheng17.h5'
    params:
        n_components = config['dim_reduction']['tsne']['n_components'],
        init = config['dim_reduction']['tsne']['init']
    output:
        ANALYSIS_OUTPUT+'/tsne/{sample}_sim_de{de_prob}_loc'+'{loc}'+'.csv'
    log:
        out = LOG_FILES+'/tsne/sample_{sample}_de{de_prob}loc_{loc}.out',
        err = LOG_FILES+'/tsne/sample_{sample}_de{de_prob}loc_{loc}.err'
    shell:
        "python scripts/apply_tsne.py {input} {output} {params.n_components} {params.init}  2> {log.err} 1> {log.out}"

rule griph:
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_de{de_prob}_loc{loc}_zheng17.h5'
    params:
        #TODO read from the config
    output:
        ANALYSIS_OUTPUT+'/griph/clusters/{sample}_sim_de{de_prob}_loc'+'{loc}'+'.csv'
    log:
        out = LOG_FILES+'/griph/sample_{sample}_de{de_prob}loc_{loc}.out',
        err = LOG_FILES+'/griph/sample_{sample}_de{de_prob}loc_{loc}.err'
    threads:8
    shell:
        "Rscript scripts/apply_griph.R --input {input} --output {output} --threads {threads} 2> {log.err} 1> {log.out}"


rule zifa:
    # by default n_blocks are estimated as n_cells/500
    input:
        SIMULATED_DATA_OUTPUT+'/{sample}_sim_de{de_prob}_loc{loc}_zheng17.h5'
    params:
        n_components = config['dim_reduction']['block_zifa']['n_components'],
        n_blocks = config['dim_reduction']['block_zifa']['n_blocks']
    output:
        ANALYSIS_OUTPUT+'/block_zifa/{sample}_sim_de{de_prob}_loc'+'{loc}'+'.csv'
    log:
        out = LOG_FILES+'/zifa/sample_{sample}_de{de_prob}loc_{loc}.out',
        err = LOG_FILES+'/zifa/sample_{sample}_de{de_prob}loc_{loc}.err'
    shell:
        "python scripts/apply_zifa.py {input} {output} {params.n_components} --n_blocks {params.n_blocks} 1> {log.out} 2> {log.err}"
