import argparse
import os.path

# parse input arguments if any are given
parser = argparse.ArgumentParser(description = 'Generate a uniform spread of amino sequences.')
parser.add_argument('--hydro',
                    help = 'the path to the table containing a map from amino acid to hydro value',
                    default = './hydro.txt')
parser.add_argument('--weight',
                    help = 'the path to the table containing a map from amino acid to weight value',
                    default = './weight.txt')
parser.add_argument('--amino',
                    help = 'the path to the table containing all possible amino acid',
                    default = './amino.txt')
parser.add_argument('--o',
                    help = 'the path to the output file')

args = parser.parse_args()


hydro_filename = args.hydro
weight_filename = args.weight
amino_filename = args.amino
out_filename = args.o

# check if all files exist and if the default output is not overwriting
# a file that already exists
