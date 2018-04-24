# sc_transcriptomics_rules.py

'''
rule cellranger_count: 
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
        out = LOG_FILES+'/{sample}/cellranger_count.out',
        err = LOG_FILES+'/{sample}/cellranger_count.err'
    shell:
        '(cd {params.cr_out}; cellranger count --id=20150306_22_27 --transcriptome={input.reference} --localcores={params.local_cores} --fastqs={input.fastqs_dir} --nosecondary 2>{log.err} 1> {log.out})'
'''

rule create_hdf5:
    input:
        #genes_file = rules.cellranger_count.output.genes_file,
        genes_file = FILTERED_GENE_BARCODE_PATH+'/filtered_gene_bc_matrices/'+T_CODE+'/genes.tsv',
        #matrix_file = rules.cellranger_count.output.matrix_file,
        matrix_file = FILTERED_GENE_BARCODE_PATH+'/filtered_gene_bc_matrices/'+T_CODE+'/matrix.mtx',
        #barcodes_file = rules.cellranger_count.output.barcodes_file
        barcodes_file = FILTERED_GENE_BARCODE_PATH+'/filtered_gene_bc_matrices/'+T_CODE+'/barcodes.tsv'
    output:
        ANALYSIS_OUTPUT+'/{sample}/raw_counts.h5'
    log:
        out = LOG_FILES+'/{sample}/create_hdf5.out',
        err = LOG_FILES+'/{sample}/create_hdf5.err'
    shell:
        'python scripts/create_hdf5.py -g {input.genes_file} -m {input.matrix_file} -b {input.barcodes_file} -o {output} 2> {log.err} 1> {log.out} '

rule filter_out_noncoding:
    input:
        rules.create_hdf5.output
    output:
        ANALYSIS_OUTPUT+'/{sample}/coding_region_only.h5'
    log:
        out = LOG_FILES+'/{sample}/filter_out_noncoding.out',
        err = LOG_FILES+'/{sample}/filter_out_noncoding/.err'
    shell:
        'Rscript scripts/select_protein_coding_genes.R --input {input} --output {output} 2> {log.err} 1> {log.out} '

rule cell_cycle_removal:
    input:
        hdf5_file = rules.filter_out_noncoding.output
    output:
        ANALYSIS_OUTPUT+'/{sample}/cell_cycle_removed.h5'
    log:
        out = LOG_FILES+'/{sample}/cell_cycle_removal.out',
        err = LOG_FILES+'/{sample}/cell_cycle_removal.err'
    shell:
        "python scripts/apply_cell_cycle_removal.py -i {input.hdf5_file} -o {output} 2> {log.err} 1> {log.out}"

rule preprocess:
    input:
        hdf5_file = rules.cell_cycle_removal.output
    params:
        n_top_genes = config['preprocess']['zheng17']['n_top_genes']
    output:
        ANALYSIS_OUTPUT+'/{sample}/zheng17.h5'
    log:
        out = LOG_FILES+'/{sample}/preprocess.out',
        err = LOG_FILES+'/{sample}/preprocess.err'
    shell:
        "python scripts/apply_preprocess.py -i {input.hdf5_file} -o {output} --n_top_genes {params.n_top_genes} 2> {log.err} 1> {log.out}"


rule phenograph:
    input:
        rules.preprocess.output
    params:
        n_neighbours = config['clustering']['phenograph']['n_neighbours'],
        log_normalize = config['clustering']['phenograph']['log_normalize']
    output:
        ANALYSIS_OUTPUT+'/{sample}/phenograph/'+'clusters.csv'
    log:
        out = LOG_FILES+'/{sample}/phenograph.out',
        err = LOG_FILES+'/{sample}/phenograph.err'
    threads: config['clustering']['phenograph']['n_jobs']
    shell:
        "python scripts/apply_phenograph.py {input} {output} {params.n_neighbours} -l {params.log_normalize} --n_threads {threads} 2> {log.err} 1> {log.out}"
