import msprime
import numpy as np

class Chromosome:
    def __init__(self, name, data, n0, n1, datatype='msprime'):

        self.name = name
        self.pos = []

        self.n0 = n0
        self.n1 = n1
        self.n2 = data.sample_size - n0 - n1

        self.H = []     # Admixed
        self.F = []     # 1st anc
        self.G = []     # 2nd anc

        if datatype == 'msprime':
            for variant in data.variants():
                self.pos.append(variant.site.position/data.sequence_length)
                self.H.append(variant.genotypes[self.n1+self.n2:])
                self.F.append(variant.genotypes[:self.n1])
                self.G.append(variant.genotypes[self.n1:self.n1+self.n2])

        self.H = np.array(self.H)
        self.F = np.array(self.F)
        self.G = np.array(self.G)


mu=1.25e-7#mutation rate per bp per generation
rho = 1.6e-9#recombination rate per bp per generation
length=int(1/rho)#1 Morgan chromosome length
chr_n = 1#number of chromosomes in each "genome"
N_diploid = 500#diploid effective population size

Tdiv = 4000#split time between two source populations

sample_size = 100#number of haploid sequences in each population


####################################################################
####################Two pulse admixture scenario####################
####################################################################

def twopulse(T1, T2, m1, m2):
    #Admixture times
    T = T1 + T2

    #Admixture proportion
    print('generating ...')
    #Three populations: populations 0 and 1 are source populations, population 2 is admixed popuation
    population_configurations = [
            msprime.PopulationConfiguration(sample_size=sample_size),
            msprime.PopulationConfiguration(sample_size=sample_size),
            msprime.PopulationConfiguration(sample_size=sample_size)]

    demographic_events=[
            msprime.SimulationModelChange(time=500, model="hudson"),
            msprime.MassMigration(T2, 2, dest=1, proportion=m1),#the more recent pulse migration from the second source population into the admixed population
            msprime.MassMigration(T, 2, dest=1, proportion=m2),#the older pulse migration from the second source population into the admixed population
            msprime.MassMigration(T, 2, dest=0, proportion=1.0),#not sure if this is a correct way to handle two MassMigrations at the same time: here we create the admixed population from the first source population
            msprime.MassMigration(Tdiv, 1, dest=0, proportion=1.0)]#split between two source populations

    data = []
    for chrom in range(chr_n):
        ts = msprime.simulate(
            Ne = N_diploid, length = length, recombination_rate=rho, mutation_rate=mu, model="dtwf", random_seed=2+69*chrom,
            population_configurations=population_configurations, demographic_events=demographic_events)
        print('chromosome {}/{} generated'.format(chrom+1, chr_n))
        data.append(Chromosome(chrom, ts, sample_size, sample_size))
        print('chromosome {}/{} submited'.format(chrom+1, chr_n))
    print('generating done!')
    return data
