import argparse
from unsupervised_methods import Zifa

parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="input matrix hdf5 file")
parser.add_argument("output_file", help="path to the output csv file")
parser.add_argument("n_components", help="the number of principal components",
        type=int)
parser.add_argument("--n_blocks", help="the number of blocks to speedup",
        type=int, default=8)
args = parser.parse_args()


def apply_zifa(input_file, output_file, n_components, n_blocks):
    
    zifa = Zifa(n_components, n_blocks)
    zifa.load_from_hdf5(input_file)
    zifa.log_normalize()
    zifa.apply()
    zifa.write_csv(output_file, dim_red_res=True)

apply_zifa(args.input_file, args.output_file, args.n_components, args.n_blocks)
