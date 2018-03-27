import argparse
from unsupervised.silhouette import Silhouette

parser = argparse.ArgumentParser()

parser.add_argument("input_file")
parser.add_argument("output_file")
parser.add_argument("h_affinity", help="Metric used to compute the linkage.")
parser.add_argument("h_linkage", help="The linkage criterion determines which distance to use between sets of observation.")
parser.add_argument("k_min", type=int)
parser.add_argument("k_max", type=int)
parser.add_argument("metric")

args = parser.parse_args()




def apply_silhouette_hierarchical(input_file, output_file, h_affinity, h_linkage, k_min, k_max, metric):

    sh = Silhouette.init_with_hierarchical(h_affinity, h_linkage, k_min, k_max, metric)
    sh.load_from_csv(input_file)
    sh.apply()
    sh.write_csv(output_file)

apply_silhouette_hierarchical(args.input_file, args.output_file, args.h_affinity, args.h_linkage, args.k_min, args.k_max, args.metric)
