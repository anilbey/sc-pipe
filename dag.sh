 snakemake --dag -s snakefile_benchmark_clust.py -j 4 | dot -Tsvg > dag.svg