import argparse
from unsupervised.simlr import Simlr

parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="input matrix hdf5 file")
parser.add_argument("output_file", help="path to the output csv file")
parser.add_argument("n_components", help="the number of latent factors for simlr",
        type=int)
parser.add_argument("pca_components", help="# PCs to use", type=int)
parser.add_argument("n_neighbours", help="#nearest neighbours", type=int)
parser.add_argument("max_iter", help="max number of iterations", type=int)
parser.add_argument("--n_threads", help="the number of threads", type=int, default=1)
args = parser.parse_args()

#  https://github.com/anilbey/SIMLR_PY simlr version is used here to support multi-threading
def apply_simlr(input_file, output_file, n_components, pca_components, n_neighbours, max_iter, threads):

    simlr = Simlr(n_components, pca_components, n_neighbours, max_iter, threads)
    simlr.load_from_hdf5(input_file)
    simlr.log_normalize()
    simlr.apply()

    simlr.write_csv(output_file)

apply_simlr(args.input_file, args.output_file, args.n_components, args.pca_components, args.n_neighbours, args.max_iter, args.n_threads)
