import msprime
from cyvcf2 import VCF
from datetime import datetime
import numpy as np

class Params:
    def __init__(self, T1, T2, M1, M2):
        self.N = 500
        self.T1 = T1
        self.T2 = T2
        self.M1 = M1
        self.M2 = M2
        self.T = self.T1 + self.T2

class Chromosome:
    def __init__(self, name, data, n0, n1, datatype='real'):

        self.name = name
        self.pos = []

        self.n0 = n0
        self.n1 = n1
        self.n2 = 0

        self.H = []     # Admixed
        self.F = []     # 1st anc
        self.G = []     # 2nd anc

        if datatype == 'msprime':
            self.n2 = data.sample_size - n0 - n1
            for variant in data.variants():
                self.pos.append(variant.site.position/data.sequence_length)
                self.H.append(variant.genotypes[self.n1+self.n2:])
                self.F.append(variant.genotypes[:self.n1])
                self.G.append(variant.genotypes[self.n1:self.n1+self.n2])

        if datatype == 'real':
            c = 0
            self.H = []
            self.F = []
            self.G = []
            self.n2 = 94 #n2 = n2

            for var in data[0]:
                self.pos.append(var.start/data[1])

                h = var.genotype.array()[0:94] #0, n
                h = (h[:,0]+h[:,1])/2

                f = var.genotype.array()[94+108:] #n, n0
                f = (f[:,0]+f[:,1])/2

                #g = var.genotype.array()[0:94] #n0 + n, n1
                #g = (g[:,0]+g[:,1])/2

                self.H.append(h)
                self.F.append(f)
                #self.G.append(g)

                c+=1
                if c%100000 == 0:
                    print(c)
            self.H = np.array(self.H)
            self.F = np.array(self.F)
            self.G = self.H
        self.H = np.array(self.H, dtype=np.byte)
        self.F = np.array(self.F, dtype=np.byte)
        self.G = np.array(self.G, dtype=np.byte)
        self.pos = np.array(self.pos, dtype=np.double)

        self.c_x = len(self.pos)


mu=1.25e-8#mutation rate per bp per generation
rho = 1.6e-9#recombination rate per bp per generation
length=int(1/rho) #1 Morgan chromosome length
N_diploid = 500#diploid effective population size

Tdiv = 4000#split time between two source populations

sample_size = 100#number of haploid sequences in each population


####################################################################
####################Two pulse admixture scenario####################
####################################################################

def twopulse(T, M, chrom):
    print('generating ...')
    #Three populations: populations 0 and 1 are source populations, population 2 is admixed popuation
    population_configurations = [
            msprime.PopulationConfiguration(sample_size=sample_size),
            msprime.PopulationConfiguration(sample_size=sample_size),
            msprime.PopulationConfiguration(sample_size=sample_size)]

    demographic_events=[
            msprime.SimulationModelChange(time=500, model="hudson"),
            msprime.MassMigration(T[1], 2, dest=1, proportion=M[1]),#the more recent pulse migration from the second source population into the admixed population
            msprime.MassMigration(T[0]+T[1], 2, dest=1, proportion=M[0]),#the older pulse migration from the second source population into the admixed population
            msprime.MassMigration(T[0]+T[1], 2, dest=0, proportion=1.0),#not sure if this is a correct way to handle two MassMigrations at the same time: here we create the admixed population from the first source population
            msprime.MassMigration(Tdiv, 1, dest=0, proportion=1.0)]#split between two source populations

    ts = msprime.simulate(
        Ne = N_diploid, length = length, recombination_rate=rho, mutation_rate=mu, model="dtwf", random_seed=2+69*chrom,
        population_configurations=population_configurations, demographic_events=demographic_events)
    chr = Chromosome(chrom, ts, sample_size, sample_size, datatype='msprime')
    print('chromosome {} submited'.format(chrom+1))
    print('generating done!')
    return chr


def twopulse_upd(T, M, chrom):
    print('generating ...')
    #Three populations: populations 0 and 1 are source populations, population 2 is admixed popuation
    population_configurations = [
            msprime.PopulationConfiguration(sample_size=sample_size),
            msprime.PopulationConfiguration(sample_size=sample_size),
            msprime.PopulationConfiguration(sample_size=sample_size)]

    dem = msprime.Demography()
    dem.add_population(name='H', initial_size=N_diploid)
    dem.add_population(name='F', initial_size=N_diploid)
    dem.add_population(name='G', initial_size=N_diploid)
    dem.add_population(name='old', initial_size=N_diploid)

    dem.add_mass_migration(time=T[1], source = "H", dest = "G", proportion=M[1])
    dem.add_admixture(time=T[0]+T[1], ancestral=["F","G"], derived="H", proportions=[1-M[0], M[0]])
    dem.add_population_split(time=Tdiv, derived=["F", "G"], ancestral="old")
    dem.sort_events()

    #print(dem)

    ts = msprime.sim_ancestry(samples={"F": sample_size/2, "G": sample_size/2, "H": sample_size/2}, demography=dem,
        sequence_length = length, recombination_rate=rho, model="dtwf", random_seed=1337+69*chrom)
    mts = msprime.sim_mutations(ts, rate=mu, random_seed = chrom*542%100+1)

    #for population in mts.populations():
    #    print(f"id={population.id}: {population.metadata}")

    chr = Chromosome(chrom, mts, sample_size, sample_size, datatype='msprime')
    print('chromosome {} submited'.format(chrom+1))
    print('generating done!')
    return chr


def realdata(file, chr_i): #this is shit
    pop_vcf = VCF(file)
    startTime = datetime.now()
    chr_data = [pop_vcf(str(chr_i + 1)), pop_vcf.seqlens[chr_i]]
    chr = Chromosome(chr_i, chr_data, 94, 108, datatype='real') #94 = n0, 108 = n1
    endTime = datetime.now()
    print('chromosome {} submited for {}'.format(chr_i+1, endTime - startTime))
    return chr
