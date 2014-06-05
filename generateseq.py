import argparse
import os.path
import sys
import csv
import numpy as np
import operator
from itertools import izip, cycle, islice


def generateDict(f, converter = float):
    ''' create a table from a csv file. The delimiters are expected to be \''\'
    '''
    table = dict()
    reader = csv.reader(f, delimiter= ',')
    for row in reader:
        table[row[0]] = converter(row[1])
    return table

def generateList(f):
    ''' create a list from a csv file. The delimiters are expected to be \''\'
    '''
    l = []
    reader = csv.reader(f, delimiter= ',')
    for row in reader:
        l.append(row[0])
    return l

def positive_int(x):
    ''' cast the input into an int and throw an argument error
        if not strictly positive
    '''
    x_int = int(x)
    if x<=0 :
        raise argparse.ArgumentTypeError("%s is an invalid strictly positive int value" % x)
    return x_int

def generateRandom(seq, size, num_samples = None):
    ''' generate a random word from a uniform distribution
    '''
    if num_samples == None:
        while True:
            yield np.random.choice(seq, size = size, replace = True)
    else:
        for i in xrange(num_samples):
            yield np.random.choice(seq, size = size, replace = True)

class BinIndex(object):
    ''' callable object that will return the index of a particular bin given
        a multi-dimensional input
    '''
    def __init__(self, minimum, maximum, nbins):
        self.min = minimum
        self.max = maximum
        self.nbins = nbins
        if len(nbins) > 1:
            self.size = reduce(operator.mul, nbins)
        else:
            self.size = nbins**min.size

    def __call__(self, x):
        norm_x = self.nbins*((x-self.min)/(self.max-self.min))
        return np.ravel_multi_index(norm_x.astype(int), mode = 'clip')

class Scorer(object):
    ''' callable object that will return the values of a particular sequence
    '''
    def __init__(self, tables):
        self.tables = tables
        self.size = len(tables)

    def __call__(self, seq):
        output = np.empty(self.size, dtype = np.double)
        for i in xrange(self.size):
            output[i] = reduce(lambda x, y: x + self.tables[i][y], seq)
        return output

# method taken from itertools recipes
def roundrobin(*iterables):
    "roundrobin('ABC', 'D', 'EF') --> A D E B F C"
    # Recipe credited to George Sakkis
    pending = len(iterables)
    nexts = cycle(iter(it).next for it in iterables)
    while pending:
        try:
            for next in nexts:
                yield next()
        except StopIteration:
            pending -= 1
            nexts = cycle(islice(nexts, pending))

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
parser.add_argument('--seqlength',
                    help = 'number of amino acid to sequence (must be non-zero, positive)',
                    default= 15,
                    type = positive_int)
parser.add_argument('--randomseed',
                    help = "set the random seed (must be integer). Use this to generate results that can be replicated. By default, a new seed will be generated using the system default behavior",
                    type = int)
parser.add_argument('--numsamples',
                    help = "set the number of random samples generated, this should be more than the number of sequences",
                    type = int,
                    default = 100000)
parser.add_argument('--maxiterations',
                    help = "set the maximum number of iterations",
                    type = int,
                    default = 1000000)
parser.add_argument('--nbins',
                    help = "set the number of bins used in the algorithm",
                    type = int,
                    default = 30)

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
    amino = np.array(generateList(f))

if hydros == None or weights == None or amino == None:
    print 'ERROR: failure to load files. Unforseen error, make sure files are properly structured.'

# set random seed if specified
if args.randomseed != None:
    np.random.seed(args.randomseed)

with open(out_filename, 'wb') as f:
    minimum = np.array([ min(hydros.values()), min(weights.values())])
    minimum *= args.seqlength

    maximum = np.array([ max(hydros.values()), max(weights.values())])
    maximum *= args.seqlength

    samples = set(generateRandom(amino, args.seqlength, args.numsamples))
    index = BinIndex(minimum, maximum, args.nbins)
    score = Scorer([hydros, weights])

    bins = { i : set() for i in xrange(index.size)}
    for s in samples:
        i = index(score(s))
        bins[i].add(s)



    counts = [[len(bins[i]), i] for i in xrange(index.size)]

    for i, new_sample in izip(xrange(args.maxiterations), generateRandom(amino, args.seqlength)):
        # remove sample sequence from largest bin
        bin_count = max(counts, key = lambda x: x[0])
        bins[bin_count[1]].remove(np.random.choice([ s for s in bins[bin_count[1]]]))
        bin_count[0] -= 1

        # add a random sequence to the appropriate bin
        bin_index = index(score(new_sample))
        bin_count = counts[bin_index]
        bins[bin_count[1]].add(new_sample)
        bin_count[0] += 1

        if i % 1000:
            bin_count = min(counts, key = lambda x: x[0])
            if bin_count[0] >= args.n/args.nbins:
                break

    results = roundrobin(*bins.values())
    for r in results:
        f.write(str(r) + '\n')










