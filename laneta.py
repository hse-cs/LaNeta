#!/usr/bin/env python3

import argparse
import sys
import time
from random import randrange
from random import randint
import numpy as np

from LaNeta.ThLd import ThLd
from LaNeta import gen

def read_simulation_setup(file):
    d = {}
    f = open(file)
    for line in f:
        (key, value) = line.split('=')
        key = key.split()[0]
        value = value.split()[0]
        d[key] = value
        print(key,'=', value)
    f.close()
    d['T']=[int(d['T1']),int(d['T2'])]
    d['M']=[float(d['M1']), float(d['M2'])]
    d['sample_sizes']=int(d['adm_size']), int(d['src1_size']), int(d['src2_size'])
    d['mu']=float(d['mu'])
    d['rho'] = float(d['rho'])
    d['N_haploid'] = float(d['N_dip'])
    d['Tdiv'] = int(d['Tdiv'])
    d['lenght_m'] = float(d['lenght'])
    d['chr_n'] = int(d['chr_n'])
    return d


parser = argparse.ArgumentParser(description='...')

#brackets
parser.add_argument('--e', '-e', nargs=1, type=float, default=0.001,
                    help='distance between brackets (default is 0.001)')
parser.add_argument('--r', '-r', nargs=1, type=float, default=0.0,
                    help='radius for bracket')

#files
parser.add_argument('--vcffile', '-vcf', nargs=1, default=None,
                    help='.vcf file that contains the admixed population')
parser.add_argument('--popfile', '-pf', nargs=1, default=None,
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
parser.add_argument("--msprime", '-ms',
                    help="Use msprime to generate data",
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

#proportions
parser.add_argument('--m1', '-m1', nargs=1, type=float, default=None,
                    help='adm. prop of the first event for the second source pop.')
parser.add_argument('--m2', '-m2', nargs=1, type=float, default=None,
                    help='adm. prop of the first event for the second source pop.')
parser.add_argument('--mt', '-mt', nargs=1, type=float, default=None,
                    help='total adm. prop. for the second source pop.')

#other
parser.add_argument('--seed', '-seed', nargs=1, type=int, default=randint(100000, 1000000),
                    help='random seed')


clargs = parser.parse_args()

#brackets
if isinstance(clargs.e, list):
    clargs.e = clargs.e[0]

if isinstance(clargs.r, list):
    clargs.r = clargs.r[0]
if clargs.r == 0:
    clargs.r = clargs.e / 2

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
if isinstance(clargs.mt, list):
    clargs.mt = clargs.mt[0]

#borders
if isinstance(clargs.max, list):
    clargs.max = clargs.max[0]
if isinstance(clargs.min, list):
    clargs.min = clargs.min[0]

#other
if isinstance(clargs.seed, list):
    clargs.seed = clargs.seed[0]


if clargs.msprime:
    ts_list = []
    print('msprime: generating...')
    d = read_simulation_setup('simulator_setup.txt')
    seed = clargs.seed
    chr_n = d['chr_n']
    print('random seed:', seed,'\n')
    for chr_i in range(chr_n):
        ts_list.append(gen.twopulse_upd(T=d['T'], M=d['M'], sample_sizes=d['sample_sizes'],
                                       mu=d['mu'], rho = d['rho'], N_haploid = d['N_haploid'],
                                       Tdiv = d['Tdiv'], lenght_m = d['lenght_m'], seed=seed*chr_n+chr_i))
        print('chrom {}/{} generated'.format(chr_i + 1, chr_n))
    exp = ThLd(data_ms=ts_list, m=clargs.mt, m1=clargs.m1, m2=clargs.m2)
else:
    exp = ThLd(popfile=clargs.popfile, vcffile=clargs.vcffile,
               pop0=clargs.pop0, pop1=clargs.pop1, pop2=clargs.pop2,
               mapfile=clargs.mapfile, m=clargs.mt, m1=clargs.m1, m2=clargs.m2)

est = exp.estimate_time(jk=clargs.jk, af=clargs.affine_term, cm_min=clargs.min,
                        cm_max=clargs.max, du=clargs.e, r=clargs.r, mt=clargs.mt)
#     exp = ThLd(popfile=clargs.popfile, vcffile=clargs.vcffile,
#                pop0=clargs.pop0, pop1=clargs.pop1, pop2=clargs.pop2,
#                mapfile=clargs.mapfile, m=clargs.mt, m1=clargs.m1, m2=clargs.m2)
#
# est = exp.estimate_time(jk=clargs.jk, af=clargs.affine_term, cm_min=clargs.min, cm_max=clargs.max,
#                         du=clargs.e, r=clargs.r, mt=clargs.mt)
print('T1, T2:', est[0], est[1])
if clargs.jk:
    print('T1: ({}, {})\nT2: ({}, {})'.format(est[2][0], est[2][1], est[3][0], est[3][1]))
