import datetime

configfile: 'config/tumor_profiler_config.json'

SAMPLE = ['melanomaS2']
T_CODE = config['transcriptome_code'] 
ANALYSIS_OUTPUT = config['analysis_output']
LOG_FILES = ANALYSIS_OUTPUT+'/log'
CELL_RANGER_OUTPUT_PATH = config['cell_ranger_output']
FILTERED_GENE_BARCODE_PATH = config['cell_ranger_filtered_matrix_path']
RUN_ID = config['unique_run_id']


include: './sc_transcriptomics_rules.py'

rule all:
    input:
        clustering_results = expand(ANALYSIS_OUTPUT+'/{sample}/phenograph/'+'clusters.csv', sample=SAMPLE)
    shell:
        'echo test rule all {input}'


