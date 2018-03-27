import argparse
from unsupervised.phenograph import Phenograph

parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="input matrix hdf5 file")
parser.add_argument("output_file", help="path to the output csv file")
parser.add_argument("n_neighbours", help="the number of neighbours",
        type=int)
parser.add_argument("--n_threads", help="the number of threads", type=int,
        default=1)
args = parser.parse_args()



def apply_pheno(input_file, output_file, n_neighbours, threads):

    pheno = Phenograph(n_neighbours, threads)
    pheno.load_from_hdf5(input_file)
    pheno.log_normalize()
    pheno.apply()

    pheno.write_csv(output_file)



apply_pheno(args.input_file, args.output_file, args.n_neighbours, args.n_threads)
