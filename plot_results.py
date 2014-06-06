import matplotlib.pyplot as plt
import os.path
import sys
import csv
import numpy as np
import operator
from itertools import izip, cycle, islice
from generateseq import generateDict, generateList, BinIndex, Scorer, generateRandom
from ast import literal_eval as make_tuple
    
hydro_filename = 'hydros.txt'
weight_filename = 'weights.txt'
amino_filename = 'codons.txt'
out_filename = 'results.txt-amino'
unfiltered_filename = 'results.txt-unfiltered'

with open(hydro_filename, 'rU') as f:
    hydros = generateDict(f)
with open(weight_filename, 'rU') as f:
    weights = generateDict(f)
with open(amino_filename, 'rU') as f:
    amino = np.array(generateList(f))

if hydros == None or weights == None or amino == None:
    print 'ERROR: failure to load files. Unforseen error, make sure files are properly structured.'

samples = []
with open(out_filename, 'rb') as f:
    for line in f:
        samples.append(make_tuple(line))
        
unfiltered_samples = []
with open(unfiltered_filename, 'rb') as f:
    for line in f:
        unfiltered_samples.append(make_tuple(line))

print 'Total number of samples: ' + str(len(samples))

score = Scorer([hydros, weights])
scores = [ tuple(score(s)) for s in samples]
plt.figure()
plt.scatter(*zip(*scores), alpha = 0.6, s = 1)

score = Scorer([hydros, weights])
scores = [ tuple(score(s)) for s in unfiltered_samples]
plt.figure()
plt.scatter(*zip(*scores), alpha = 0.6, s = 2)
