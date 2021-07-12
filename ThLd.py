#!/usr/bin/env python3

import argparse
import sys
import time
from random import randrange
import numpy as np

from clThLd import ThLd
import gen

parser = argparse.ArgumentParser(description='...')


parser.add_argument('--p', '-p', nargs=1, type=int, default=300,
                    help='number of discrete bins (default is 300)')
parser.add_argument('--chrn', '-cn', nargs=1, type=int, default=20,
                    help='number of chromosome pairs (default is 23)')
parser.add_argument('--sampleSizeH', '-Hs', nargs=1, type=int, default=100,
                    help='number of sample (default is 100)')
parser.add_argument('--sampleSizeF', '-Fs', nargs=1, type=int, default=100,
                    help='number of sample (default is 100)')
parser.add_argument('--sampleSizeG', '-Gs', nargs=1, type=int, default=100,
                    help='number of sample (default is 100)')
parser.add_argument('--seed', '-seed', nargs=1, type=int, default=None,
                    help='random seed')
parser.add_argument('--filename', '-f', nargs=1, default=None,
                    help='vcf file with data')
parser.add_argument("--msprime", '-ms',
                    help="Use msprime to generate data",
                    action="store_true")


clargs = parser.parse_args()

if isinstance(clargs.p, list):
    clargs.p = clargs.p[0]
if isinstance(clargs.chrn, list):
    clargs.chrn = clargs.chrn[0]
if isinstance(clargs.sampleSizeH, list):
    clargs.sampleSizeH = clargs.sampleSizeH[0]
if isinstance(clargs.sampleSizeF, list):
    clargs.sampleSizeF = clargs.sampleSizeF[0]
if isinstance(clargs.sampleSizeG, list):
    clargs.sampleSizeG = clargs.sampleSizeG[0]
if isinstance(clargs.seed, list):
    clargs.seed = clargs.seed[0]
if isinstance(clargs.filename, list):
    clargs.filename = clargs.filename[0]


chroms = []
if clargs.msprime:
    print('msprime: generating...')
    params = np.loadtxt('tap.txt')
    TrueTs = params[0]
    TrueMs = params[1]
    for chr_i in range(clargs.chrn):
        print('\nchrom {}/{}'.format(chr_i + 1, clargs.chrn))
        #self.chrom = gen.twopulse(self.TrueTs, self.TrueMs, chr_i)
        chroms.append(gen.twopulse_upd(TrueTs, TrueMs, chr_i))
else:
    for chr_i in range(clargs.chrn):
        print('real data')
        print('\nchrom {}/{}'.format(chr_i + 1, clargs.chrn))
        chroms.append(gen.realdata(clargs.filename, chr_i))
exp = ThLd(clargs.p, chroms)
exp.alg()
exp.estimate()
