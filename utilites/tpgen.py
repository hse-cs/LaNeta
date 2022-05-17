#
# You can use this python3 script for
# generating two pulse model data
# in apropriate format.
#
# Dependencies:
#   msprime
#   vcftools
#   tabix
#   bgzip
#   merger.sh
#
# sintaxis:
# python3 tpgen.py $sample $chr_n
# sample - integer, sample, used as uniq random seed and output directory
# chr_n - integer, number of chromosomes in sample.
#

import msprime

import subprocess
import sys
import os

# parameters:
t2 = 10
t1 = 5
m1 = 0.2
m2 = 0.2
N = 1000
sample_size=100
seed = 213142332
#

sample = int(sys.argv[1])
chr_n = int(sys.argv[2])

print('generating two pulse model data:')
print('sample:', sample)
print('chr. num.:', chr_n)

def twopulse_upd(T=[10,10], M=[0.2, 0.2], sample_sizes=[100, 100, 100], mu=1.25e-8, rho = 1.6e-9, N_haploid = [1000, 1000, 1000, 1000], Tdiv = 4000, lenght_m = 1, seed=1):
    length=int(lenght_m/rho)

    dem = msprime.Demography()
    dem.add_population(name='H', initial_size=N_haploid[0])
    dem.add_population(name='F', initial_size=N_haploid[1])#, default_sampling_time=T[0]+T[1])
    dem.add_population(name='G', initial_size=N_haploid[2])
    dem.add_population(name='old', initial_size=N_haploid[3])

    dem.add_mass_migration(time=T[1] + 1, source = "H", dest = "G", proportion=M[1])
    dem.add_admixture(time=T[0]+T[1] + 1, ancestral=["F","G"], derived="H", proportions=[1-M[0], M[0]])
    dem.add_population_split(time=Tdiv, derived=["F", "G"], ancestral="old")
    dem.sort_events()

    ts = msprime.sim_ancestry(
        samples={"H": sample_sizes[0], "F": sample_sizes[1], "G": sample_sizes[2]},
        demography=dem, sequence_length = length, recombination_rate=rho, ploidy=2,
        model=[msprime.DiscreteTimeWrightFisher(duration=50), msprime.StandardCoalescent(duration=3950)],
        #model='hudson',
        random_seed=seed)
    mts = msprime.sim_mutations(ts, rate=mu, random_seed=seed)
    return mts

admixed_names = ['H_' + str(i) for i in range(0, sample_size)]
reference_1_names = ['F_' + str(i) for i in range(0, sample_size)]
reference_2_names = ['G_' + str(i) for i in range(0, sample_size)]
ind_names = admixed_names + reference_1_names + reference_2_names

try:
    os.mkdir(str(sample))
except OSError as error:
    print('sample already exists, rewriting...')
    subprocess.call(['rm', '-r', str(sample)])
    os.mkdir(str(sample))

with open(str(sample)+'/sample_'+str(sample)+".map", "w") as mapfile:
    mapfile.close()

with open(str(sample)+'/sample_'+str(sample)+".pop", "w") as popfile:
    for name in admixed_names:
        ent = name + ' H\n'
        popfile.write(ent)
    for name in reference_1_names:
        ent = name + ' F\n'
        popfile.write(ent)
    for name in reference_2_names:
        ent = name + ' G\n'
        popfile.write(ent)
    popfile.close()

for chr_i in range(1, chr_n + 1):
    ts = twopulse_upd(T=[t1,t2], M=[m1, m2],
             sample_sizes=[sample_size, sample_size, sample_size],
             mu=1.25e-8, rho = 1.6e-9,
             N_haploid = [N, N, N, N],
             Tdiv = 4000, lenght_m = 1, seed=seed + sample*chr_n + chr_i)
    with open(str(sample)+'/sample_'+str(sample)+'_chr_'+str(chr_i)+".vcf", "wt") as vcf_file:
        ts.write_vcf(vcf_file, individual_names = ind_names, contig_id=chr_i)
        vcf_file.close()

    pos = ts.tables.sites.position

    with open(str(sample)+'/sample_'+str(sample)+".map", "a") as mapfile:
        for x in pos:
            ent = str(chr_i) + ' . ' + str(x) + ' ' + str(x/ts.sequence_length*100) + '\n'
            mapfile.write(ent)
        mapfile.close()

    print(f'chr: {chr_i}/{chr_n} - done!')

print('merging chromosomes in one .vcf ...')
subprocess.call(["./merger.sh", str(sample)])
print('Done!')
