import argparse
import logging
from unsupervised.phenograph import Phenograph
from utils import str2bool

logging.basicConfig(level=logging.DEBUG)


parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="input matrix hdf5 file")
parser.add_argument("output_file", help="path to the output csv file")
parser.add_argument("n_neighbours", help="the number of neighbours",
        type=int)
parser.add_argument("--n_threads", help="the number of threads", type=int,
        default=1)
parser.add_argument("-l", "--log_normalize", help="Boolean switch for the log normalization", type=str2bool)
args = parser.parse_args()



def apply_pheno(input_file, output_file, n_neighbours, threads, log_normalize):

    logging.debug('log_normalization value: ' + str(log_normalize))
    logging.debug('type of log_normalization value: ' + str(type(log_normalize)))

    pheno = Phenograph(n_neighbours, threads)
    pheno.load_from_hdf5(input_file)
    if log_normalize:
        pheno.log_normalize()
    else:
        pass
    pheno.apply()

    pheno.write_csv(output_file)



apply_pheno(args.input_file, args.output_file, args.n_neighbours, args.n_threads, args.log_normalize)
