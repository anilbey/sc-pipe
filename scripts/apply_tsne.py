import argparse
from unsupervised_methods import Tsne

parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="input matrix hdf5 file")
parser.add_argument("output_file", help="path to the output csv file")
parser.add_argument("n_components", help="the number of components",
        type=int)
parser.add_argument("init")
parser.add_argument("--n_threads", help="the number of threads", type=int,
        default=1)
args = parser.parse_args()



def apply_tsne(input_file, output_file, n_components, init, threads):

    tsne = Tsne(n_components, init)
    tsne.load_from_hdf5(input_file)
    tsne.log_normalize()
    tsne.apply()
    tsne.write_csv(output_file, dim_red_res=True)

apply_tsne(args.input_file, args.output_file, args.n_components, args.init, args.n_threads)
