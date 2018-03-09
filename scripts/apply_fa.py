import argparse
from unsupervised_methods import FactorAnalysis

parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="input matrix hdf5 file")
parser.add_argument("output_file", help="path to the output csv file")
parser.add_argument("n_components", help="the number of latent factors",
        type=int)
parser.add_argument("--n_threads", help="the number of threads", type=int,
        default=1)
args = parser.parse_args()



def apply_fa(input_file, output_file, n_components, threads):

    fa = FactorAnalysis(n_components)
    fa.load_from_hdf5(input_file)
    fa.log_normalize()
    fa.apply()

    fa.write_csv(output_file, dim_red_res=True)

apply_fa(args.input_file, args.output_file, args.n_components, args.n_threads)
