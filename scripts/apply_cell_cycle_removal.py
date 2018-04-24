import argparse
from cell_cycle import remove_cell_cycle

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", required=True, help="path to the input hdf5 file")
parser.add_argument("-o", "--output_file", required=True, help="path to the output hdf5 file")
args = parser.parse_args()

remove_cell_cycle(args.input_file, args.output_file)
