import argparse
from preprocess import preprocess

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", required=True, help="path to the input hdf5 file")
parser.add_argument("-o", "--output_file", required=True, help="path to the output hdf5 file")
parser.add_argument("--n_top_genes", required=True, type=int, help="the number of genes to be selected")
args = parser.parse_args()

preprocess(args.input_file, args.output_file, args.n_top_genes)
