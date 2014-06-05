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

def generateRandom(seq, size, num_samples = None, tester = lambda x: True):
    ''' generate a random word from a uniform distribution. Tester returns true
        if a given sequence is legal. This generator will never terminate if
        it is impossible to satisfy tester.
    '''
    if num_samples == None:
        while True:
            new_sample = tuple(np.random.choice(seq, size = size, replace = True))
            while not tester(new_sample):
                new_sample = tuple(np.random.choice(seq, size = size, replace = True))
            yield new_sample
    else:
        for i in xrange(num_samples):
            new_sample = tuple(np.random.choice(seq, size = size, replace = True))
            while not tester(new_sample):
                new_sample = tuple(np.random.choice(seq, size = size, replace = True))
            yield new_sample

def filtersequences(samples, scorer, tester):
    ''' filter a list of sequences. tester compares to sequences are returns
        true if they are not too close. The returned list is a legal set of
        sequence according to the test. There are no optimality guarantees.
    '''
    sample_score = sorted([ (s, scorer(s)) for s in samples], key = lambda x: x[1][0])
    sample_set = [sample_score[0]]
    for i in xrange(1, len(sample_score)):
        if np.all([tester(s[1], sample_score[i][1]) for s in sample_set]):
            sample_set.append(sample_score[i])
    return zip(*sample_set)[0]

class BinIndex(object):
    ''' callable object that will return the index of a particular bin given
        a multi-dimensional input
    '''
    def __init__(self, minimum, maximum, nbins):
        self.min = minimum
        self.max = maximum
        self.nbins = nbins
        if isinstance(nbins, int) or len(nbins) == 1:
            self.size = nbins**self.min.size
            if not isinstance(nbins, int):
                self.nbins = [nbins[0]]*self.min.size
            else:
                self.nbins = [nbins]*self.min.size
        else:
            self.size = reduce(operator.mul, nbins)

    def __call__(self, x):
        norm_x = self.nbins*((x-self.min)/(self.max-self.min))
        return np.ravel_multi_index(norm_x.astype(int),
                                    dims = self.nbins,
                                    mode = 'clip')

class Scorer(object):
    ''' callable object that will return the values of a particular sequence
    '''
    def __init__(self, tables):
        self.tables = tables
        self.size = len(tables)

    def __call__(self, seq):
        output = np.empty(self.size, dtype = np.double)
        for i in xrange(self.size):
            output[i] = sum(map(lambda y: self.tables[i][y], seq))
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

if __name__ == '__main__':
    # parse input arguments if any are given
    parser = argparse.ArgumentParser(description = 'Generate a uniform spread of amino sequences.')
    parser.add_argument('--hydros',
                        help = 'the path to the table containing a map from amino acid to hydro value',
                        default = './hydros.txt')
    parser.add_argument('--weights',
                        help = 'the path to the table containing a map from amino acid to weight value',
                        default = './weights.txt')
    parser.add_argument('--codons',
                        help = 'the path to the table containing all possible amino acid',
                        default = './codons.txt')
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
                        type = positive_int,
                        default = 100000)
    parser.add_argument('--maxiterations',
                        help = "set the maximum number of iterations",
                        type = positive_int,
                        default = 100000)
    parser.add_argument('--nbins',
                        help = "set the number of bins used in the algorithm",
                        type = positive_int,
                        metavar = 'N',
                        nargs = '+',
                        default = [30])

    args = parser.parse_args()


    hydro_filename = args.hydros
    weight_filename = args.weights
    codon_filename = args.codons
    out_filename = args.o

    # check if all files exist and if output exits, check if overwriting
    # is allowed
    filename_error = False
    error_msg = 'ERROR: File \'{}\' not found for {}. If you haven\'t explicitly specified a filepath a default file path was used'
    overwrite_error_msg = 'ERROR: File \'{}\' exists and will not be overwritten unless the \'--f\' argument is given'
    if not os.path.exists(hydro_filename):
        filename_error = True
        print error_msg.format(hydro_filename, 'hydro values')

    if not os.path.exists(weight_filename):
        filename_error = True
        print error_msg.format(weight_filename, 'weight values')

    if not os.path.exists(codon_filename):
        filename_error = True
        print error_msg.format(codon_filename, 'codon values')

    if os.path.exists(out_filename) and not args.f:
        filename_error = True
        print overwrite_error_msg.format( out_filename)

    if filename_error:
        sys.exit(1)

    # load tables and list
    hydros = None
    weights = None
    amino = None
    codons = None

    with open(hydro_filename, 'rU') as f:
        hydros = generateDict(f)
    with open(weight_filename, 'rU') as f:
        weights = generateDict(f)
    with open(codon_filename, 'rU') as f:
        codons = generateDict(f, str)
    with open(codon_filename, 'rU') as f:
        amino = np.array(generateList(f))

    if hydros == None or weights == None or amino == None:
        print 'ERROR: failure to load files. Unforseen error, make sure files are properly structured.'

    # set random seed if specified
    if args.randomseed != None:
        np.random.seed(args.randomseed)

    with open(out_filename, 'wb') as f, open(out_filename+'-amino', 'wb') as fa:
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

            index_to_remove = np.random.choice(len(bins[bin_count[1]]))
            items = list(bins[bin_count[1]])

            bins[bin_count[1]].remove(items[index_to_remove])
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
            f.write(str(tuple([ codons[x] for x in r])) + '\n')
            fa.write(str(r) + '\n')










