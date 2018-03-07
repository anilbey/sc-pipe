import datetime

configfile: 'config/config.json'

# global variables
SAMPLE = ['melanomaS2']
T_CODE = config['transcriptome_code'] # e.g. hg19
ANALYSIS_OUTPUT = config['real_data_output']
LOG_FILES = ANALYSIS_OUTPUT+'/log'
CELL_RANGER_OUTPUT_PATH = config['cell_ranger_output']
HDF5_OUTPUT = 'hdf5_data'
RUN_ID = config['unique_run_id']


# TODO: have a method to generate the unique_run_id by using the sample name, date and time info
#def u_run_id():
#    return datetime.datetime.now().strftime ("%Y%m%d_%H_%M")

rule all:
    input:
        clustering_results = expand(ANALYSIS_OUTPUT+'/{sample}/phenograph/'+'clusters.csv', sample=SAMPLE)
#TODO input function        fastqs = get_all_fastqs(config['input_fastqs'])
    shell:
        'echo test rule all {input}'
'''
rule cellranger_count: # (parallel)
    input:
        fastqs_dir = config['input_fastqs']+'/{sample}',
        reference = config['reference_transcriptome']
    params:
        cr_out = CELL_RANGER_OUTPUT_PATH+'/{sample}',
        id = RUN_ID,
        local_cores = config['cellranger_count']['local_cores']
    output:
        genes_file = CELL_RANGER_OUTPUT_PATH + '/{sample}/' +RUN_ID+ '/outs/filtered_gene_bc_matrices/'+T_CODE+'/genes.tsv',
        matrix_file = CELL_RANGER_OUTPUT_PATH + '/{sample}/' +RUN_ID+'/outs/filtered_gene_bc_matrices/'+T_CODE+'/matrix.mtx',
        barcodes_file = CELL_RANGER_OUTPUT_PATH + '/{sample}/' +RUN_ID+'/outs/filtered_gene_bc_matrices/'+T_CODE+'/barcodes.tsv'
    log:
        out = LOG_FILES+'/cellranger_count/sample_{sample}.out',
        err = LOG_FILES+'/cellranger_count/sample_{sample}.err'
    shell:
        '(cd {params.cr_out}; cellranger count --id=20150306_22_27 --transcriptome={input.reference} --localcores={params.local_cores} --fastqs={input.fastqs_dir} --nosecondary 2>{log.err} 1> {log.out})'
'''
rule create_hdf5:
    input:
        #genes_file = rules.cellranger_count.output.genes_file,
        genes_file = CELL_RANGER_OUTPUT_PATH+'/20150306_22_27/outs/filtered_gene_bc_matrices/'+T_CODE+'/genes.tsv',
        #matrix_file = rules.cellranger_count.output.matrix_file,
        matrix_file = CELL_RANGER_OUTPUT_PATH+'/20150306_22_27/outs/filtered_gene_bc_matrices/'+T_CODE+'/matrix.mtx',
        #barcodes_file = rules.cellranger_count.output.barcodes_file
        barcodes_file = CELL_RANGER_OUTPUT_PATH+'/20150306_22_27/outs/filtered_gene_bc_matrices/'+T_CODE+'/barcodes.tsv'
    output:
        HDF5_OUTPUT+'/raw_{sample}.h5'
    log:
        out = LOG_FILES+'/create_hdf5/sample_{sample}.out',
        err = LOG_FILES+'/create_hdf5/sample_{sample}.err'
    shell:
        'python scripts/create_hdf5.py -g {input.genes_file} -m {input.matrix_file} -b {input.barcodes_file} -o {output} 2> {log.err} 1> {log.out} '

rule preprocess_zheng17:
    input:
        hdf5_file = rules.create_hdf5.output
    params:
        n_top_genes = config['preprocess']['zheng17']['n_top_genes']
    output:
        HDF5_OUTPUT+'/zheng17_{sample}.h5'
    log:
        out = LOG_FILES+'/preprocess_zheng17/sample_{sample}.out',
        err = LOG_FILES+'/preprocess_zheng17/sample_{sample}.err'
    shell:
        "python scripts/preprocess_zheng17.py -i {input.hdf5_file} -o {output} --n_top_genes {params.n_top_genes} 2> {log.err} 1> {log.out}"

rule phenograph:
    input:
        #rules.preprocess_zheng17.output
        HDF5_OUTPUT+'/zheng17_{sample}.h5'
    params:
        n_neighbours = config['clustering']['phenograph']['n_neighbours']
    output:
        ANALYSIS_OUTPUT+'/{sample}/phenograph/'+'clusters.csv'
    log:
        out = LOG_FILES+'/phenograph/sample_{sample}.out',
        err = LOG_FILES+'/phenograph/sample_{sample}.err'
    threads: config['clustering']['phenograph']['n_jobs']
    shell:
        "python scripts/apply_phenograph.py {input} {output} {params.n_neighbours} --n_threads {threads} 2> {log.err} 1> {log.out}"
