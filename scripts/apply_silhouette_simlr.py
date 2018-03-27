import argparse
from unsupervised.silhouette import Silhouette

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)

parser.add_argument("--n_components", help="the number of latent factors for simlr",
        type=int, required=True)
parser.add_argument("--pca_components", help="# PCs to use", type=int, required=True)
parser.add_argument("--n_neighbours", help="#nearest neighbours", type=int, required=True)
parser.add_argument("--max_iter", help="max number of iterations", type=int, required=True)
parser.add_argument("--n_threads", help="the number of threads", type=int, default=1)
parser.add_argument("--k_min", type=int, required=True)
parser.add_argument("--k_max", type=int, required=True)
parser.add_argument("--metric", required=True)

args = parser.parse_args()

def apply_silhouette_simlr(input_file, output_file, n_components, pca_components, n_neighbours, max_iter, k_min, k_max, metric, threads):

    ss = Silhouette.init_with_simlr(n_components, pca_components, n_neighbours, max_iter, threads, k_min, k_max, metric)
    ss.load_from_hdf5(input_file)
    ss.log_normalize()
    ss.apply()
    ss.write_csv(output_file)


apply_silhouette_simlr(args.input, args.output, args.n_components, args.pca_components, args.n_neighbours, args.max_iter, args.k_min, args.k_max, args.metric, args.n_threads)
