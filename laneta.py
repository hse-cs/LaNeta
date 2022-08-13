#!/usr/bin/env python3

import argparse

from LaNeta.ThLd import LaNeta


parser = argparse.ArgumentParser(description='...')

#brackets
parser.add_argument('--binsize', '-b', nargs=1, type=float, default=0.001,
                    help='distance between brackets (default is 0.001)')
parser.add_argument('--binradius', '-r', nargs=1, type=float, default=0,
                    help='radius for bracket, default is binsize/2')

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

parser.add_argument('--chromosomes', '-c', default=[], nargs='+',
                    help='name of chromosomes to proceed (default: all contigs from .vcf header are used)')

#flags
parser.add_argument("--haplotype", '-ht',
                    help="Use if you use haplotype data instead of genotype",
                    action="store_true")
parser.add_argument('--jk', '-jk',
                    help="Estimate confidence intervals with jackknife",
                    action="store_true")
parser.add_argument('--nmt', '-nmt',
                    help="To use mt only for allele frequencies.",
                    action="store_true")

#borders
parser.add_argument('--max', '-max', nargs=1, type=float, default=30,
                    help='max d and ds for estimation in cantimograns')
parser.add_argument('--min', '-min', nargs=1, type=float, default=0.5,
                    help='min d and ds for estimation in cantimograns')

#parameters
parser.add_argument('--mt', '-mt', nargs=1, type=float, default=None,
                    help='total adm. prop. for the second source pop.')
parser.add_argument('--m1', '-m1', nargs=1, type=float, default=None,
                    help='adm. prop of the first event for the second source pop.')
parser.add_argument('--m2', '-m2', nargs=1, type=float, default=None,
                    help='adm. prop of the first event for the second source pop.')
parser.add_argument('--t1', '-t1', nargs=1, type=float, default=None,
                    help='adm. time of the first event.')
parser.add_argument('--t2', '-t2', nargs=1, type=float, default=None,
                    help='adm. time of the first event.')


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

#flags
if isinstance(clargs.jk, list):
    clargs.jk = clargs.jk[0]
if isinstance(clargs.nmt, list):
    clargs.nmt = clargs.nmt[0]
if isinstance(clargs.haplotype, list):
    clargs.haplotype = clargs.haplotype[0]

print('Starting LaNeta...')

# защита
if clargs.binsize < 0 or clargs.binradius < 0:
    print('Wrong bin format!')
elif clargs.pop0 == None or (clargs.pop1 == None and clargs.pop2 == None):
    print('Admixed population and at least one source population are required!')
elif (clargs.pop1 == None or clargs.pop2 == None) and clargs.mt == None:
    print('If only one source population is avalable total admixture proportion (-mt) is requared!')
elif clargs.mapfile == None or clargs.popfile == None or clargs.vcffile == None:
    print('Please provide all requared input files! (-vcf, -m, -p)')
else:
    exp = LaNeta(popfile=clargs.popfile, vcffile=clargs.vcffile,
                 pop0=clargs.pop0, pop1=clargs.pop1, pop2=clargs.pop2,
                 mapfile=clargs.mapfile, gt=not clargs.haplotype, contig_names=clargs.chromosomes)

    parameters = exp.estimateTwoPulse(cm_max=clargs.max,
                                      bin_size=clargs.binsize, bin_radius=clargs.binradius,
                                      cm_min=clargs.min, jk=clargs.jk,
                                      M1=clargs.m1, M2=clargs.m2,
                                      T1=clargs.t1, T2=clargs.t2,
                                      Mt=clargs.mt, nmt=clargs.nmt)
    print()
    exp.print_parameters(jk=clargs.jk)
