import argparse
from unsupervised.pca import Pca

parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="input matrix hdf5 file")
parser.add_argument("output_file", help="path to the output csv file")
parser.add_argument("n_components", help="the number of principal components",
        type=int, default=2)
parser.add_argument("--n_threads", help="the number of threads", type=int,
        default=1)
args = parser.parse_args()


def apply_pca(input_file, output_file, n_components, threads):

    pca = Pca(n_components)
    pca.load_from_hdf5(input_file)
    pca.log_normalize()
    pca.apply()

    pca.write_csv(output_file, dim_red_res=True)

apply_pca(args.input_file, args.output_file, args.n_components, args.n_threads)
