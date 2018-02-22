import argparse
from unsupervised_methods import Silhouette

parser = argparse.ArgumentParser()

parser.add_argument("input_file")
parser.add_argument("output_file")
parser.add_argument("n_init", type=int)
parser.add_argument("n_jobs", type=int)
parser.add_argument("k_min", type=int)
parser.add_argument("k_max", type=int)
parser.add_argument("metric")

args = parser.parse_args()




def apply_silhouette_kmeans(input_file, output_file, n_init, n_jobs, k_min, k_max, metric):

    skm = Silhouette.init_with_kmeans(n_init, n_jobs, k_min, k_max, metric)
    skm.load_from_csv(input_file)
    skm.apply()
    skm.write_csv(output_file)

apply_silhouette_kmeans(args.input_file, args.output_file, args.n_init, args.n_jobs, args.k_min, args.k_max, args.metric)
