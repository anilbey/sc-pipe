library(griph)
library(rhdf5)
library(argparse)

parser <- ArgumentParser(description= 'simulates data')
parser$add_argument('--input', required=TRUE)
parser$add_argument('--output', required=TRUE)
parser$add_argument('--threads', required=TRUE, type='integer')
args = parser$parse_args()



apply_griph = function(input_file, output_file, threads)
{
    h5f = H5Fopen(input_file)
    matrix = h5f$matrix
    H5Fclose(h5f)
    matrix = log1p(matrix)
    res = griph_cluster(matrix, ncores=threads, plot_ = FALSE)

    write.table(x=res$MEMB,file=output_file, sep=',', row.names=FALSE, col.names=FALSE)

}

apply_griph(args$input, args$output, args$threads)
