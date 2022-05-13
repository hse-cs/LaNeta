#!/usr/bin/env python3

import argparse
import sys
import time
from random import randrange
from random import randint
import numpy as np

from LaNeta.ThLd import LaNeta


parser = argparse.ArgumentParser(description='...')

#brackets
parser.add_argument('--binsize', '-b', nargs=1, type=float, default=0.001,
                    help='distance between brackets (default is 0.001)')
parser.add_argument('--binradius', '-r', nargs=1, type=float, default=0.0,
                    help='radius for bracket')

#files
parser.add_argument('--vcffile', '-vcf', nargs=1, default=None,
                    help='.vcf file that contains the admixed population')
parser.add_argument('--popfile', '-p', nargs=1, default=None,
                    help='.txt file that contains samples with indicated population')
parser.add_argument("--mapfile", '-m',
                    help=".txt file with chromosome name(1-22), var id, var position(bp), var position(cm)")
parser.add_argument('--pop0', '-p0', nargs=1, default=None,
                    help='name of the admixed population in popfile')
parser.add_argument('--pop1', '-p1', nargs=1, default=None,
                    help='name of the first source population in popfile')
parser.add_argument('--pop2', '-p2', nargs=1, default=None,
                    help='name of the second source in popfile')

#flags
parser.add_argument("--genotype", '-gt',
                    help="Use if you use genotype data instead of haplotype",
                    action="store_true")
parser.add_argument("--affine_term", '-af',
                    help="Estimate affine term that due to population substructure",
                    action="store_true")
parser.add_argument('--jk', '-jk',
                    help="Estimate confidence intervals with jackknife",
                    action="store_true")

#borders
parser.add_argument('--max', '-max', nargs=1, type=int, default=20,
                    help='max d and ds for estimation in cantimograns')
parser.add_argument('--min', '-min', nargs=1, type=int, default=1,
                    help='min d and ds for estimation in cantimograns')

#parameters
parser.add_argument('--m1', '-m1', nargs=1, type=float, default=None,
                    help='adm. prop of the first event for the second source pop.')
parser.add_argument('--m2', '-m2', nargs=1, type=float, default=None,
                    help='adm. prop of the first event for the second source pop.')
parser.add_argument('--t1', '-t1', nargs=1, type=float, default=None,
                    help='adm. time of the first event.')
parser.add_argument('--t2', '-t2', nargs=1, type=float, default=None,
                    help='adm. time of the first event.')
parser.add_argument('--mt', '-mt', nargs=1, type=float, default=None,
                    help='total adm. prop. for the second source pop.')


clargs = parser.parse_args()

#brackets
if isinstance(clargs.binsize, list):
    clargs.binsize = clargs.binsize[0]

if isinstance(clargs.binradius, list):
    clargs.binradius = clargs.binradius[0]
if clargs.binradius == 0:
    clargs.binradius = clargs.binsize / 2

#files
if isinstance(clargs.vcffile, list):
    clargs.vcffile = clargs.vcffile[0]
if isinstance(clargs.mapfile, list):
    clargs.mapfile = clargs.mapfile[0]
if isinstance(clargs.popfile, list):
    clargs.popfile = clargs.popfile[0]
if isinstance(clargs.pop0, list):
    clargs.pop0 = clargs.pop0[0]
if isinstance(clargs.pop1, list):
    clargs.pop1 = clargs.pop1[0]
if isinstance(clargs.pop2, list):
    clargs.pop2 = clargs.pop2[0]

#proportions
if isinstance(clargs.m1, list):
    clargs.m1 = clargs.m1[0]
if isinstance(clargs.m2, list):
    clargs.m2 = clargs.m2[0]
if isinstance(clargs.t1, list):
    clargs.t1 = clargs.t1[0]
if isinstance(clargs.t2, list):
    clargs.t2 = clargs.t2[0]
if isinstance(clargs.mt, list):
    clargs.mt = clargs.mt[0]

#borders
if isinstance(clargs.max, list):
    clargs.max = clargs.max[0]
if isinstance(clargs.min, list):
    clargs.min = clargs.min[0]

print('Starting LaNeta...')
exp = LaNeta(popfile=clargs.popfile, vcffile=clargs.vcffile,
             pop0=clargs.pop0, pop1=clargs.pop1, pop2=clargs.pop2,
             mapfile=clargs.mapfile, gt=clargs.genotype)

parameters = exp.estimateTwoPulse(cm_max=20, freq_filter=0,
                                  bin_size=0.001, bin_radius=0.0005,
                                  cm_min=0.5,
                                  M1=clargs.m1, M2=clargs.m2,
                                  T1=clargs.t1, T2=clargs.t2,
                                  Mt=clargs.mt)


print('[T1, T2, M1, M2] =', parameters)
if clargs.jk:
    print('T1: ({}, {})\nT2: ({}, {})'.format(est[2][0], est[2][1], est[3][0], est[3][1]))
