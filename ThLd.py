#!/usr/bin/env python3

import argparse
import sys
import time
from random import randrange
from random import randint
import numpy as np

from clThLd import ThLd
import gen

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


parser.add_argument('--e', '-e', nargs=1, type=float, default=0.001,
                    help='bracket width (default is 0.001)')
parser.add_argument('--admixed_vcf', '-a', nargs=1, default=None,
                    help='.vcf file that contains the admixed population')
parser.add_argument('--source1_vcf', '-s1', nargs=1, default=None,
                    help='.vcf file that contains the first source population')
parser.add_argument('--source2_vcf', '-s2', nargs=1, default=None,
                    help='.vcf file that contains the second source population')
parser.add_argument("--mapfile", '-m',
                    help=".txt file with chromosome name(1-22), var id, var position(bp), var position(cm)")
parser.add_argument("--msprime", '-ms',
                    help="Use msprime to generate data",
                    action="store_true")
parser.add_argument("--affine_term", '-af',
                    help="Estimate affine term (for real data) that due to population substructure",
                    action="store_true")
parser.add_argument('--seed', '-seed', nargs=1, type=int, default=randint(100000, 1000000),
                    help='random seed')
parser.add_argument('--cm', '-cm', nargs=1, type=int, default=20,
                    help='max d and ds for estimation in cantimograns')
parser.add_argument('--m1', '-m1', nargs=1, type=float, default=None,
                    help='adm. prop of the first event for the second source pop.')
parser.add_argument('--m2', '-m2', nargs=1, type=float, default=None,
                    help='adm. prop of the first event for the second source pop.')
parser.add_argument('--mt', '-mt', nargs=1, type=float, default=None,
                    help='total adm. prop. for the second source pop.')

clargs = parser.parse_args()
if isinstance(clargs.e, list):
    clargs.e = clargs.e[0]
if isinstance(clargs.seed, list):
    clargs.seed = clargs.seed[0]
if isinstance(clargs.admixed_vcf, list):
    clargs.admixed_vcf = clargs.admixed_vcf[0]
if isinstance(clargs.source1_vcf, list):
    clargs.source1_vcf = clargs.source1_vcf[0]
if isinstance(clargs.source2_vcf, list):
    clargs.source2_vcf = clargs.source2_vcf[0]
if isinstance(clargs.mapfile, list):
    clargs.mapfile = clargs.mapfile[0]
if isinstance(clargs.m1, list):
    clargs.m1 = clargs.m1[0]
if isinstance(clargs.m2, list):
    clargs.m2 = clargs.m2[0]
if isinstance(clargs.mt, list):
    clargs.mt = clargs.mt[0]
if isinstance(clargs.cm, list):
    clargs.cm = clargs.cm[0]


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
    exp = ThLd(data_ms=ts_list)
else:
    exp = ThLd(vcf0=clargs.admixed_vcf, vcf1=clargs.source1_vcf, vcf2=clargs.source2_vcf, mapfile=clargs.mapfile)
    if clargs.affine_term:
        exp.start_affine(max_c_af=10000)

exp.start(du = clargs.e)

if clargs.m1!=None and clargs.m2!=None:
    exp.M = [clargs.m1, clargs.m2]
    print(exp.M, clargs.cm)
    x = exp.estimate_ls(cm=clargs.cm)
elif clargs.Mt != None:
    exp.Mt = clargs.mt
    x = exp.estimate_ls_one_prop(cm=clargs.cm)

print(x, exp.M)
