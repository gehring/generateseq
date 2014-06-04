import argparse
import os.path
import sys
import csv
import numpy as np
from itertools import izip

def generateDict(file):
    table = dict()
    reader = csv.reader(file, delimiter= ',')
    for row in reader:
        table[row[0]] = row[1]
    return table

def generateList(file):
    l = []
    reader = csv.reader(file, delimiter= ',')
    for row in reader:
        l.append(row[0])
    return l

def positive_int(x):
    x_int = int(x)
    if x<=0 :
        raise argparse.ArgumentTypeError("%s is an invalid strictly positive int value" % x)
    return x_int

def generateRandom(seq, length):
    while True:
        yield np.random.choice(seq, size = length, replace = True)

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
                    help = 'the path to the output file',
                    default = './results.txt')
parser.add_argument('--f',
                    help = 'allow the output file to be overwritten',
                    action="store_true")
parser.add_argument('--n',
                    help = 'number of sequences to generate (must be non-zero, positive)',
                    default=1000,
                    type = positive_int)
parser.add_argument('--seq-length',
                    help = 'number of amino acid to sequence (must be non-zero, positive)',
                    default= 15,
                    type = positive_int)
parser.add_argument('--random-seed',
                    help = "set the random seed (must be integer). Use this to generate results that can be replicated. By default, a new seed will be generated using the system default behavior",
                    type = int)

args = parser.parse_args()


hydro_filename = args.hydro
weight_filename = args.weight
amino_filename = args.amino
out_filename = args.o

# check if all files exist and if output exits, check if overwriting
# is allowed
filename_error = False
error_msg = 'ERROR: File \'%s\' not found for %s. If you haven\'t explicitly specified a filepath a default file path was used'
overwrite_error_msg = 'ERROR: File \'%s\' exists and will not be overwritten unless the \'--f\' argument is given'
if not os.path.exists(hydro_filename):
    filename_error = True
    print error_msg % hydro_filename, 'hydro values'

if not os.path.exists(weight_filename):
    filename_error = True
    print error_msg % weight_filename, 'weight values'

if not os.path.exists(amino_filename):
    filename_error = True
    print error_msg % amino_filename, 'amino values'

if os.path.exists(out_filename) and not args.f:
    filename_error = True
    print overwrite_error_msg % out_filename

if filename_error:
    sys.exit(1)

# load tables and list
hydros = None
weights = None
amino = None

with open(hydro_filename, 'rb') as f:
    hydros = generateDict(f)
with open(weight_filename, 'rb') as f:
    weights = generateDict(f)
with open(amino_filename, 'rb') as f:
    amino = generateList(f)

if hydros == None or weights == None or amino == None:
    print 'ERROR: failure to load files. Unforseen error, make sure files are properly structured.'

# set random seed if secified
if args.random-seed != None:
    np.random.seed(args.random-seed)




