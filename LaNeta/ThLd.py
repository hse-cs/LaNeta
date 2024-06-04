from LaNeta.fastc import FFT, thld_6, binning_deltas, estimate_frequencies_c, estimate_deltas_c, \
noFFT, binning_bins_int, estimate_affine_term, estimate_coef_slow, estimate_affine_term_v2
from LaNeta.fastc import read_vcf_genotype, read_positions
import numpy as np
from datetime import datetime
from scipy.optimize import minimize
from scipy.optimize import least_squares

from cyvcf2 import VCF


def get_vcf_seqlen(vcffile, name):
    vcf = VCF(vcffile)
    sequence_length = 0
    for var in vcf(name):
        sequence_length += 1
    vcf.close()
    return sequence_length

def read_pos(mapfile, chr_name):
    position_bp, position, sequence_length = read_positions(mapfile, chr_name)
    pos = np.array(position)
    pos = np.sort(pos)
    pos_bp = np.array(position_bp)
    pos_bp = np.sort(pos_bp)
    return pos_bp, pos, sequence_length

def smart_add(a, b):
    if a.shape[0]>b.shape[0]:
        c = np.zeros((a.shape[0], a.shape[0]))
        c += a
        c[:b.shape[0], :b.shape[0]] += b
        return c
    else:
        c = np.zeros((b.shape[0], b.shape[0]))
        c += b
        c[:a.shape[0], :a.shape[0]] += a
        return c

def mkarr(ls):
    shape = 0
    for arr in ls:
        if arr.shape[0] > shape:
            shape = arr.shape[0]
    Arr = np.zeros((len(ls), shape, shape))
    for i in range(len(ls)):
        Arr[i, :ls[i].shape[0], :ls[i].shape[0]] = ls[i]
    return Arr

def get_individuals_nodes(ts, population_id, sampled=False):
    nodes = [] # nodes of samples
    for node in ts.nodes(): # they are contained in return of ts.nodes()
        if node.population == population_id and node.individual != -1:
            nodes.append(node.id)
    return nodes


def read_msprime(data_input, chr_i, c_max = None, pos_read=True, gt=False):
    data = Chromosome()

    ts_list = data_input[0]
    pop0, pop1, pop2 = data_input[1], data_input[2], data_input[3]

    ts = ts_list[chr_i]
    data.n0 = len(ts.samples(pop0))# * 10
    data.n1 = len(ts.samples(pop1))
    data.n2 = len(ts.samples(pop2))

    adm_samples_id = get_individuals_nodes(ts, pop0)
    src1_samples_id = get_individuals_nodes(ts, pop1)
    src2_samples_id = get_individuals_nodes(ts, pop2)

    pos = np.empty(ts.sites().length, dtype=np.double)
    H = np.empty((ts.sites().length, data.n0), dtype=np.byte)     # Admixed
    F = np.empty((ts.sites().length, data.n1), dtype=np.byte)     # 1st anc
    G = np.empty((ts.sites().length, data.n2), dtype=np.byte)    # 2nd anc
    #print(H.shape, F.shape, G.shape)
    c = ts.sites().length
    if data.n1 == 0:
        H = np.empty((ts.sites().length, data.n0 - data.n0//2), dtype=np.byte)     # 1st anc
        F = np.empty((ts.sites().length, data.n0//2), dtype=np.byte)    # 2nd anc
        data.missing_ref = 1
    if data.n2 == 0:
        H = np.empty((ts.sites().length, data.n0 - data.n0//2), dtype=np.byte)     # 1st anc
        G = np.empty((ts.sites().length, data.n0//2), dtype=np.byte)    # 2nd anc
        data.missing_ref = 2
    if pos_read:
        pos = ts.tables.sites.position / ts.sequence_length
        gt_matrix = ts.genotype_matrix()
        if data.n1 == 0:
            F = gt_matrix[:, adm_samples_id[:data.n0//2]]
            H = gt_matrix[:, adm_samples_id[data.n0//2:]]
            G = gt_matrix[:, src2_samples_id]
        elif data.n2 == 0:
            G = gt_matrix[:, adm_samples_id[:data.n0//2]]
            H = gt_matrix[:, adm_samples_id[data.n0//2:]]
            F = gt_matrix[:, src1_samples_id]
        else:
            G = gt_matrix[:, src2_samples_id]
            F = gt_matrix[:, src1_samples_id]
            H = gt_matrix[:, adm_samples_id]
    else:
        c_total = ts.get_num_sites()
        print('max:', c_max, 'total:', c_total)
        rng = np.random.default_rng()
        inds = np.sort(rng.choice(c_total, size=c_max, replace=False))

        id = 0
        j = 0
        H = H[:c_max]
        F = F[:c_max]
        G = G[:c_max]
        for variant in ts.variants():
            if j != len(inds) and id == inds[j]:
                if data.n1 == 0:
                    F[j] = variant.genotypes[: data.n0 // 2]
                    H[j] = variant.genotypes[data.n0 // 2: data.n0]
                    G[j] = variant.genotypes[data.n0:]
                elif data.n2 == 0:
                    G[j] = variant.genotypes[: data.n0 // 2]
                    H[j] = variant.genotypes[data.n0 // 2: data.n0]
                    F[j] = variant.genotypes[data.n0:]
                else:
                    G[j] = variant.genotypes[data.n0 + data.n1 : data.n0 + data.n1 + data.n2]
                    F[j] = variant.genotypes[data.n0 : data.n0 + data.n1]
                    H[j] = variant.genotypes[: data.n0]
                j += 1
            id += 1
        c = c_max


    data.n0 = H.shape[1]
    data.n1 = F.shape[1]
    data.n2 = G.shape[1]

    # # make gt data:
    if gt:
        for i in range(0, data.n0, 2):
            H[:, i//2] = (H[:, i] + H[:, i+1])

        for i in range(0, data.n1, 2):
            F[:, i//2] = (F[:, i] + F[:, i+1])

        for i in range(0, data.n2, 2):
            G[:, i//2] = (G[:, i] + G[:, i+1])

        data.n0//=2
        data.n1//=2
        data.n2//=2
    # #
    # print(data.n0, data.n1, data.n2, H.dtype)
    data.G = np.array(G[:, :data.n2], dtype=np.byte, order='C')
    data.H = np.array(H[:, :data.n0], dtype=np.byte, order='C')
    data.F = np.array(F[:, :data.n1], dtype=np.byte, order='C')
    data.pos = pos[:c]
    data.sequence_length = c
    #print('pa',H.shape)

    # if gt:
    #     data.n0*=2
    #     data.n1*=2
    #     data.n2*=2
    data.isEmpty = False
    return data

def read_pop(vcffile, popfile, populations):

    vcf = VCF(vcffile)
    samples = vcf.samples
    vcf.close()

    H_pop, F_pop, G_pop = populations
    H_samples = []
    F_samples = []
    G_samples = []

    f = open(popfile)
    for line in f:
        sample, pop = line.split()
        sample = sample#+'_'+sample #AAAAAAAAAA ANASTASIA IMPROVE IT PLZ!
        if pop == H_pop:
            if sample in samples:
                H_samples.append(samples.index(sample))
        if pop == F_pop:
            if sample in samples:
                F_samples.append(samples.index(sample))
        if pop == G_pop:
            if sample in samples:
                G_samples.append(samples.index(sample))

    return np.array(H_samples, dtype=np.int), np.array(F_samples, dtype=np.int), np.array(G_samples, dtype=np.int)

def read_real(files, chr_name, c_max = None, pos_read=True, **kwargs):
    #open vcf
    data = Chromosome()
    vcffile = files[0]
    popfile = files[1]
    mapfile = files[2]
    populations = files[3], files[4], files[5]


    H_samples, F_samples, G_samples = read_pop(vcffile, popfile, populations)
    data.n1 = len(F_samples)
    data.n2 = len(G_samples)
    data.n0 = len(H_samples)

    #read vcf and txt
    F, G, H, pos = None, None, None, None
    if data.n1 == 0:
        F_samples = H_samples[:data.n0//2]
        H_samples = H_samples[data.n0//2:]
        data.n1 = data.n0//2
        data.n0 = data.n0 - data.n1
        data.missing_ref = 1
    if data.n2 == 0:
        G_samples = H_samples[:data.n0//2]
        H_samples = H_samples[data.n0//2:]
        data.n2 = data.n0//2
        data.n0 = data.n0 - data.n2
        data.missing_ref = 2

    vcf = VCF(vcffile)
    if pos_read:
        print(f'CHR: {chr_name}', end=' ')
        pos_bp, pos, pos_entry_count = read_pos(mapfile, chr_name) #chr_i+1 because .txt contains names of chomosomes
        if pos_entry_count == 0:
            print('- Empty .map!')
            data.isEmpty = True
            return data
        vcf_entry_count = get_vcf_seqlen(vcffile ,chr_name)
        if vcf_entry_count == 0:
            print('- Empty!')
            data.isEmpty = True
            return data
        print()
        # print(f'length: {pos[pos_entry_count-1]-pos[0]} Morgans.')
        print(f'variants count (.vcf): {vcf_entry_count}')
        if pos_entry_count != vcf_entry_count:
            print(f'variants count (.pos): {pos_entry_count}')
            print(f'Number of lines for CHR {chr_name} in .vcf and .pos is different!')
            if pos_entry_count < vcf_entry_count:
                print('Using only variants from .pos file...')
            else:
                print('.pos is not valid!')
        sequence_length = min(pos_entry_count, vcf_entry_count)
        H, F, G = read_vcf_genotype(vcf(chr_name), pos_bp, sequence_length, H_samples, F_samples, G_samples)
        data.isEmpty = False
    else:
        _, _, c_total = read_pos(mapfile, c_max, chr_name)
        print('max:', c_max, 'total:', c_total)
        rng = np.random.default_rng()
        inds = np.sort(rng.choice(c_total, size=c_max, replace=False))
        H, F, G = read_vcf_genotype(vcf(chr_name), c_max, H_samples, F_samples, G_samples, inds=inds)
        c = c_max
    vcf.close()

    data.G = G
    data.H = H
    data.F = F
    #pos = np.linspace(0, 2.86279, c)
    #print(np.array(pos))
    data.pos = pos
    data.sequence_length = sequence_length
    return data


class Chromosome:
    '''
    Represents chromosome data

    H, F, G - source populations 0, 1, 2.

    n0, n1, n2 - sample sizes of populations 0, 1, 2.

    missing_ref - missing source population, if 0 every source population
                  is available.
    '''
    def __init__(self): # data = [H_gen, F_gen, G_gen, seqlens]
        self.H, self.F, self.G = None, None, None
        self.H_freq, self.F_freq, self.G_freq = None, None, None
        self.n0, self.n1, self.n2 = 0, 0, 0
        self.missing_ref = 0
        self.isEmpty = True

    def estimate_frequencies(self):
        self.H_freq, self.F_freq, self.G_freq = estimate_frequencies_c(
                                                self.H, self.n0,
                                                self.F, self.n1,
                                                self.G, self.n2,
                                                self.sequence_length
                                                )

    def swap_ref(self):
        "Swap admixed and 1 available source population."
        if self.missing_ref == 2:
            self.H, self.G = self.G, self.H
            self.n0, self.n2 = self.n2, self.n0
            self.H_freq, self.G_freq = self.G_freq, self.H_freq
        elif self.missing_ref == 1:
            self.H, self.F = self.F, self.H
            self.n0, self.n1 = self.n1, self.n0
            self.H_freq, self.F_freq = self.F_freq, self.H_freq
#        print('ref and adm swaped!')


class Weighted_LD_Algorithm:
    def __init__(self, WLD):
        self.WLD = WLD

    def estimate_coef_naive(self):
        WLD = self.WLD
        eps = WLD.eps
        self.numerator, self.denumerator = estimate_coef_slow(WLD.data.H,   WLD.data.n0,
                                              WLD.data.F,   WLD.data.n1,
                                              WLD.data.G,   WLD.data.n2,
                                              WLD.data.pos, WLD.data.sequence_length,
                                              WLD.bin_size,
                                              eps)

        self.numerator, self.denumerator = np.array(self.numerator), np.array(self.denumerator)
        WLD.denumerator = self.denumerator
        WLD.numerator = self.numerator
        return self.get_coef()


    def estimate_affine_term(self):
        WLD = self.WLD
        H_chroms = np.zeros((len(WLD.chroms), WLD.max_c_af, WLD.chroms[0].n0), dtype=np.double)
        F_chroms = np.zeros((len(WLD.chroms), WLD.max_c_af, WLD.chroms[0].n1), dtype=np.double)
        G_chroms = np.zeros((len(WLD.chroms), WLD.max_c_af, WLD.chroms[0].n2), dtype=np.double)
        for i in range(len(WLD.chroms)):
            H_chroms[i, :, :] = WLD.chroms[i].H
            F_chroms[i, :, :] = WLD.chroms[i].F
            G_chroms[i, :, :] = WLD.chroms[i].G

        at = estimate_affine_term_v2(H_chroms, F_chroms, G_chroms, len(WLD.chroms), WLD.max_c_af,
                                  WLD.chroms[0].n0, WLD.chroms[0].n1, WLD.chroms[0].n2,
                                  WLD.m, WLD.chroms[0].missing_ref)

        return at

    def estimate_bins_int(self):
        WLD = self.WLD
        # if WLD.data.missing_ref != 0:
        #     print('missing freq:', WLD.data.missing_ref)
        self.WLD.b, self.WLD.c = binning_bins_int(
                                   WLD.data.H, WLD.data.n0,
                                   self.WLD.Delta_freq,
                                   WLD.data.H_freq,
                                   WLD.data.F_freq,
                                   WLD.data.G_freq,
                                   WLD.freq_filter,
                                   WLD.data.pos, WLD.data.sequence_length,
                                   WLD.bin_size, WLD.bin_radius)
        i = -1
        while self.WLD.c[i] == 0:
            i-=1
        self.border = len(self.WLD.c) + i + 1
        #WLD.b_int.append(b)
        return self.WLD.b, self.WLD.c

    def estimate_deltas(self):
        WLD = self.WLD
        deltas, self.WLD.Delta_freq = estimate_deltas_c(
                              WLD.data.H_freq,
                              WLD.data.F_freq,
                              WLD.data.G_freq,
                              WLD.freq_filter, WLD.data.sequence_length,
                              m=WLD.Mt, missing_ref=WLD.data.missing_ref)
        return deltas

    def estimate_deltas_numer(self): #dont use
        print('estimating weights (FFT).. ')
        WLD = self.WLD
        deltas_bin = np.array(binning_deltas(WLD.data.H,   WLD.data.n0,
                                             WLD.data.F,   WLD.data.n1,
                                             WLD.data.G,   WLD.data.n2,
                                             WLD.data.pos, WLD.data.sequence_length,
                                             WLD.bin_size,
                                             m=WLD.m, missing_ref=WLD.data.missing_ref) )
        deltas_numer = FFT(deltas_bin[:self.border])[1:self.border, 1:self.border]
        return deltas_numer

    def estimate_numerator(self):
        WLD = self.WLD
        b = self.WLD.b #no epsilon during binning
        #b = b / (WLD.data.n0*WLD.data.n1*WLD.data.n2) #back to freq
        #print(b[:,:self.border].shape)
        #b[:, self.c<100] = 0 # FOR REAL
        numerator = FFT(b[:,:self.border])[1:self.border, 1:self.border]
        #self.numerator = np.round(self.numerator) / ( (WLD.data.n0*WLD.data.n1*WLD.data.n2)**3 )
        #WLD.numerator = self.numerator
        return numerator

    def estimate_numerator_slow(self):
        WLD = self.WLD
        b = self.estimate_bin_b_int()
        #print(b[:,:self.border].shape)
        self.numerator = self.numerator / ( (WLD.data.n0*WLD.data.n1*WLD.data.n2)**3 )
        self.numerator = noFFT(b[:,:self.border], b.shape[0], b.shape[1])
        return self.numerator

    def estimate_denumerator_slow(self):
        WLD = self.WLD
        c = self.estimate_bin_c()
        c = c[:self.border].reshape(1, -1)
        c = np.array(c, dtype=np.double)
        self.numerator = noFFT(c, c.shape[0], c.shape[1])
        return self.numerator

    def estimate_denumerator(self):
        WLD = self.WLD
        c = self.WLD.c
        #print(c.shape)
        denumerator = FFT(c[:self.border])[1:self.border, 1:self.border] #zeros in the end of array
        denumerator = np.round(denumerator) #kill machine epsilon (everewhere int)
        return denumerator

    def get_delta_fft(self):
        WLD = self.WLD
        delta = np.divide(WLD.delta_numer, WLD.denumerator, out=np.zeros_like(WLD.delta_numer), where=WLD.denumerator!=0)
        return delta

    def get_coef(self):
        WLD = self.WLD
        n = WLD.data.n0
        coef = np.divide(WLD.numerator, WLD.denumerator, out=np.zeros_like(WLD.numerator), where=WLD.denumerator!=0)*(n/((n-1)*(n-2))) #default out is zeros, where divisor has zero, there is no devision
        if WLD.gt:
            coef *= 4
        return coef

    def estimate_weighted_ld_coef(self):
        print('estimating coef.. ', end='')
        WLD = self.WLD
        def direct():
            self.denumerator = self.estimate_denumerator()
            self.numerator = self.estimate_numerator()
            coef = self.get_coef()
            return coef
        coef = direct()
        if WLD.unb:
            WLD.data.swap_ref()
            coef += direct()
            coef /= 2
        print('done!')
        return coef


class Weighted_LD:
    def __init__(self,
                 data_ms=None, mapfile=None, vcffile=None, # data_ms - list of ts
                 popfile=None, pop0=None, pop1=None, pop2=None,
                 gt=False, contig_names=None):

        self.calculated = False

        self.ld_filename = 'ld.txt'
        if vcffile != None:
            self.ld_filename = vcffile.split('.')[0] + '.ld'

        self.coef = None
        self.delta = None
        self.ld = None

        self.at = 0

        #data load init
        self.one_ref=False
        self.gt = gt
        if data_ms==None:
            vcf = VCF(vcffile)
            if contig_names == None:
                self.chr_names = vcf.seqnames # 20
                print('The number of chromosomes (contigs in .vcf HEADER):', end = ' ')
            else:
                if len(set(contig_names).intersection(set(vcf.seqnames))) == 0:
                    print('None of your contig names are found in .vcf header! Be careful!')
                self.chr_names = contig_names
                print('The number of chromosomes (contigs in provided list):', end=' ')
            self.chr_n = len(self.chr_names)
            print(self.chr_n)
            if not gt:
                print('Haplotypes data provided!')
            vcf.close()

            if pop1 == None or pop2 == None:
                    self.one_ref = True
            self.read = read_real
            self.data_input = [vcffile, popfile, mapfile,
                               pop0, pop1, pop2]
            pop_0_samples, pop_1_samples, pop_2_samples = read_pop(self.data_input[0], self.data_input[1], self.data_input[3:])
            if self.one_ref:
                print('Only one source population is available')
            print('\nTwo Pulse Model:')
            if pop0 != None:
                print(pop0, f'- admixed population! Sample size: {len(pop_0_samples)}')
            else:
                print('Admixed population is missing! Please provide it!')
            to_print_pop1 = pop1 if pop1!= None else 'UNKNOWN'
            print(to_print_pop1, f'- admixed one time! Sample size: {len(pop_1_samples)}')
            to_print_pop2 = pop2 if pop2!= None else 'UNKNOWN'
            print(to_print_pop2, f'- admixed two times! Sample size: {len(pop_2_samples)}')
            if pop1 == None and pop2 == None or len(pop_1_samples) == 0 and len(pop_2_samples) == 0:
                print('Both source populations are missing! Please, provide at least one!')
        else:
            self.chr_n = len(data_ms)
            self.chr_names = range(self.chr_n)
            if len(data_ms[0].samples(1)) == 0 or len(data_ms[0].samples(2)) == 0:
                self.one_ref = True
            self.read = read_msprime
            self.data_input = [data_ms, pop0, pop1, pop2]




    def estimate(self, bin_size=0.001, bin_radius=0.0005, freq_filter=10,
                 at=False, Mt=None):
        self.Mt = Mt if Mt != None else -1
        self.freq_filter = freq_filter/100
        self.coef = np.empty(0)
        #self.delta_numer = np.empty(0)
        self.bin_size = bin_size
        self.bin_radius = bin_radius


        g_num = 0 #chr loop

        self.numerators = []
        self.denumerators = []
        self.deltas = []

        chr_n = self.chr_n
        #self.delta_numers = []

        if self.one_ref and self.Mt <= 0:
            print("Estimating WLD error!")
            print("Only one source pop. available and there is no total adm. proportion!")
        print('\nestimating weighted LD...')
        for chr_name in self.chr_names:
            t0 = datetime.now()
            self.data = self.read(self.data_input, chr_name, gt=self.gt)
            if self.data.isEmpty:
                chr_n -= 1
                continue
            print(f'{datetime.now() - t0} - data read.')
            t0 = datetime.now()

            self.data.estimate_frequencies()
            alg = Weighted_LD_Algorithm(self)
            delta_i = alg.estimate_deltas()

            alg.estimate_bins_int()
            denumerator_i = alg.estimate_denumerator()
            numerator_i = alg.estimate_numerator()

            #delta_numer_i = alg.estimate_deltas_numer()

            if self.one_ref:
                self.data.swap_ref()

                delta_i += alg.estimate_deltas()
                delta_i /= 2
                alg.estimate_bins_int()
                numerator_i += alg.estimate_numerator()
                numerator_i /= 2

                #delta_numer_i = smart_add(delta_numer_i, alg.estimate_deltas_numer())
                #delta_numer_i /= 2

            #print('{}/{}'.format(g_num, chr_n), datetime.now()-t0, 'delta_numer_mean:', delta_numer_i.mean())
            #self.delta_numer = smart_add(self.delta_numer, alg.estimate_deltas_numer())
            print(f'{datetime.now() - t0} - estimation.')
            #t0 = datetime.now()
            self.numerators.append(numerator_i)
            self.denumerators.append(denumerator_i)
            #self.delta_numers.append(delta_numer_i)
            self.deltas.append(delta_i)
            g_num+=1
            #print(f'{datetime.now() - t0} - enter.')
            print('{}/{}'.format(g_num, chr_n))#, 'delta_i:', delta_i)
            #break

        if self.numerators != [] and self.denumerators != [] and self.deltas != []:
            self.numerators = mkarr(self.numerators)
            self.denumerators = mkarr(self.denumerators)
            self.deltas = np.array(self.deltas)
            #self.delta_numers = mkarr(self.delta_numers)
            #final estimation
            alg = Weighted_LD_Algorithm(self)
            #print(g_num)

            self.numerator = self.numerators.sum(axis=0)
            self.denumerator = self.denumerators.sum(axis=0)

            self.numerator /= g_num
            self.denumerator /= g_num

            #self.delta_numer /= g_num
            #self.delta_fft = alg.get_delta_fft()
            #print(self.delta_fft)

            self.coef = alg.get_coef()

            self.delta = self.deltas.sum()
            self.delta /= g_num

            if at:
                self.at = self.estimate_affine_term()
                self.coef -= self.at

            self.ld = self.coef/self.delta**3
            np.savetxt(self.ld_filename, self.ld)
            self.calculated = True
        else:
            print('No WLD calculated!')

    def estimate_affine_term(self):
        pass

    def jackknife(self, chr_i):
        alg = Weighted_LD_Algorithm(self)
        self.numerator = np.delete(self.numerators, chr_i, axis=0).sum(axis=0)
        self.denumerator = np.delete(self.denumerators, chr_i, axis=0).sum(axis=0)
        self.numerator /= self.chr_n - 1
        self.denumerator /= self.chr_n - 1
        self.coef = alg.get_coef()
        if self.at != 0:
            self.coef -= self.at
        # self.delta_numer = np.delete(self.delta_numers, chr_i, axis=0).sum(axis=0)
        # self.delta_numer /= self.chr_n - 1
        # self.delta_fft = alg.get_delta_fft()
        self.delta = np.delete(self.deltas, chr_i, axis=0).sum()
        self.delta /= self.chr_n - 1
        self.ld = self.coef/self.delta**3

    def load_weighted_ld(self, path, cm_scale=1):
        self.coef = np.load(path)
        self.bin_size = cm_scale

    def load_delta(self, delta):
        self.delta = delta

    def load_ld(self, path, cm_scale=1):
        self.ld = np.load(path)
        self.bin_size = cm_scale

    def getld(self, cm_min=0.5, cm_max=20):
        if self.delta == None:
            print('No ld estimated')
            return -1
        P_ = int(cm_max/(100*self.bin_size))
        _P = int(cm_min/(100*self.bin_size))
        #delta = self.delta_fft[_P:P_, _P:P_]
        return self.ld[_P:P_, _P:P_]

    def getwld(self, cm_min=0.5, cm_max=20):
        if self.delta == None:
            print('No wld estimated')
            return -1
        P_ = int(cm_max/(100*self.bin_size))
        _P = int(cm_min/(100*self.bin_size))
        #delta = self.delta_fft[_P:P_, _P:P_]
        return self.coef[_P:P_, _P:P_]

    def getdelta(self):
        if self.delta == None:
            print('No delta estimated')
            return -1
        return self.delta


class LaNeta:
    """
    User interface class
    """
    def __init__(self, data_ms=None, mapfile=None, #data_ms - list of ts
                 popfile=None, pop0=None, pop1=None,
                 pop2=None, vcffile=None, gt=False, contig_names=None):

        self.cm_min = 0.5
        self.cm_max = 30

        self.jk_done = False

        self.Ms = []
        self.parameters = None
        if contig_names == []:
            contig_names = None
        self.WLD = Weighted_LD(data_ms=data_ms, mapfile=mapfile, vcffile=vcffile,
                               popfile=popfile, pop0=pop0, pop1=pop1, pop2=pop2,
                               gt=gt, contig_names=contig_names)

    def debug(self):
        self.data = Chromosome()
        one_ref=False
        if self.data_ms==None:
            vcf = VCF(self.vcffile)
            chr_n = len(vcf.seqlens) # 20
            vcf.close()
            if self.pop1 == None or self.pop2 == None:
                one_ref = True
            read = read_real
            raw = [self.vcffile, self.popfile, self.mapfile, self.pop0, self.pop1, self.pop2]

        #check vcf


        #check pos
        self.chr_n = chr_n
        g_num = 0 #chr loop

        vcf = VCF(self.vcffile)
        chr_names = vcf.seqnames
        vcf.close()

        for chr_name in chr_names:
            vcf_len = get_vcf_seqlen(self.vcffile, chr_name)
            print('vcf len:',vcf_len)
            pos_bp, pos, txt_len = read_pos(self.mapfile, vcf_len*2, chr_name)
            print('txt len:',txt_len)
            if vcf_len != txt_len:
                print('number of variants in .txt and .vcf is different!')



        #check pop

        #check calculations

    def getwld(self, cm_min=None, cm_max=None):
        if cm_min == None:
            cm_min = self.cm_min
        if cm_max == None:
            cm_max = self.cm_max
        return self.WLD.getwld(cm_min=cm_min, cm_max=cm_max)

    def getld(self, cm_min=None, cm_max=None):
        if cm_min == None:
            cm_min = self.cm_min
        if cm_max == None:
            cm_max = self.cm_max
        return self.WLD.getld(cm_min=cm_min, cm_max=cm_max)

    def load_ld(self, path, cm_scale=1):
        self.WLD.load_ld(path, cm_scale=1)
    # def start_slow(self):
    #         g_num = 0
    #         for chr_i in self.chroms:
    #             alg = Algorithm(self)
    #             self.data = chr_i
    #             self.denumerator += alg.estimate_denumerator_slow()
    #             self.numerator += alg.estimate_numerator_slow()
    #             self.deltas += alg.estimate_deltas()
    #             g_num+=1
    #             print('{}/{}'.format(g_num, len(self.chroms)))
    #         alg = Algorithm(self)
    #         alg.numerator = self.numerator
    #         alg.denumerator = self.denumerator
    #         self.coef = alg.get_coef()/2
    #         self.deltas /= g_num
    #
    # def start_slow_real(self, bin_size=0.001, eps=0.000005):
    #     self.eps = eps
    #     self.bin_size = bin_size
    #     if self.file!=None and self.mapfile!=None:
    #
    #         vcf = VCF(self.file)
    #         mapm = open(self.mapfile, 'r')
    #         g_num = 0
    #
    #         print('sequence: '+vcf.seqnames[0])
    #         self.chroms = []
    #         self.chroms = [gen.realdata(vcf, mapm, 0, self.user_n)]
    #         self.chr_n = len(self.chroms)
    #         self.data = self.chroms[0]
    #         g_num+=1
    #
    #         alg = Algorithm(self)
    #         self.coef = alg.estimate_coef_naive()
    #         self.delta += alg.estimate_deltas()
    #         print('{}/{}'.format(g_num, len(vcf.seqnames)))
    #
    #         vcf.close()
    #         mapm.close()


    # def start_affine(self, max_c_af=10000):
    #     self.max_c_af = max_c_af
    #     g_num = 0
    #     self.chroms = []
    #
    #     if self.data_ms==None:
    #         vcf = VCF(self.vcffile)
    #         chr_names = vcf.seqnames # 20
    #         chr_n = len(chr_names)
    #         vcf.close()
    #         if self.pop1 == None or self.pop2 == None:
    #             one_ref = True
    #         read = read_real
    #         raw = [self.vcffile, self.popfile, self.mapfile, self.pop0, self.pop1, self.pop2]
    #     else:
    #         chr_n =len(self.data_ms)
    #         read = read_msprime
    #         raw = self.data_ms
    #
    #     for chr_name in chr_names:
    #         data = read(raw, chr_name, pos_read=False, c_max = max_c_af)
    #         self.chroms.append(data)
    #         g_num+=1
    #         if g_num == 22:
    #             break
    #
    #     alg = Algorithm(self)
    #     self.at = alg.estimate_affine_term()

    def estimateTwoPulse(self, wld_estimation=True, jk=False, af=False,
                               M1=None, M2=None, T1=None, T2=None,
                               bin_size=0.001, bin_radius=0.0005,
                               cm_min=None, cm_max=None, at=False, Mt=None,
                               freq_filter=-1, nmt=False):
        if cm_min != None:
            self.cm_min = cm_min
        else:
            cm_min = self.cm_min
        if cm_max != None:
            self.cm_max = cm_max
        else:
            cm_max = self.cm_max

        self.WLD.estimate(bin_size=bin_size, bin_radius=bin_radius,
                          freq_filter=freq_filter, Mt=Mt, at=at)
        if self.WLD.calculated:
            print('\nestimating parameters...')
            print(f'cm min: {cm_min}\ncm max: {cm_max}')
            parameters = self.estimate_parameters(T1=T1, T2=T2, M1=M1, M2=M2,
                                                  Mt=Mt, nmt=nmt,
                                                  cm_scale=bin_size,
                                                  cm_min=cm_min, cm_max=cm_max,
                                                  save=True)
            self.parameters = parameters
            if jk:
                if self.WLD.chr_n < 2:
                    print('Less that 2 chromosomes provided, can\'t run jackknife!')
                else:
                    parameters_do = np.zeros((self.WLD.chr_n, len(parameters)))
                    self._parameters = np.zeros(len(parameters))
                    self.parameters_ = np.zeros(len(parameters))
                    print('\nrunning jackknife...')
                    for chr_i in range(self.WLD.chr_n):
                        print(f'estimating without CHR {self.WLD.chr_names[chr_i]}...')
                        self.WLD.jackknife(chr_i)
                        parameters_do[chr_i] = self.estimate_parameters(T1=T1, T2=T2, M1=M1, M2=M2,
                                                                        Mt=Mt, silent=True,
                                                                        cm_scale=bin_size, nmt=nmt,
                                                                        cm_min=cm_min, cm_max=cm_max,
                                                                        rand_mean=parameters)
                    #print(parameters_do)
                    parameters_ps = self.WLD.chr_n * np.array(parameters) - (self.WLD.chr_n-1) * parameters_do
                    # _parameters = parameters_do.mean(axis=0) - 1.96*np.sqrt((1/self.WLD.chr_n)*parameters_do.var(axis=0))
                    # parameters_ = parameters_do.mean(axis=0) + 1.96*np.sqrt((1/self.WLD.chr_n)*parameters_do.var(axis=0))
                    self.jk_bias = -((self.WLD.chr_n-1) * np.array(parameters) - (self.WLD.chr_n-1) * parameters_do.mean(axis=0))
                    self.jk_var = parameters_ps.var(axis=0)/self.WLD.chr_n
                    self.jk_interval_wide = 1.96*np.sqrt((1/self.WLD.chr_n)*self.jk_var)
                    self.jk_done = True
                    return self.parameters, self.jk_bias, self.jk_var, self.jk_interval_wide
            return self.parameters
    def print_parameters(self, jk=False):
        if self.parameters != None:
            print('t1\tt2\tm1\tm2')
            print(f'{self.parameters[0]}\t{self.parameters[1]}\t{self.parameters[2]}\t{self.parameters[3]}\t estimation')
            if jk and self.jk_done:
                print(f'{self.jk_bias[0]}\t{self.jk_bias[1]}\t{self.jk_bias[2]}\t{self.jk_bias[3]}\tjk bias')
                print(f'{self.jk_var[0]}\t{self.jk_var[1]}\t{self.jk_var[2]}\t{self.jk_var[3]}\tjk var')
                print(f'{self.jk_interval_wide[0]}\t{self.jk_interval_wide[1]}\t{self.jk_interval_wide[2]}\t{self.jk_interval_wide[3]}\tjk 95% conf. interval (+/-)')
        else:
            print('No parameters estimated!')

    def estimate_parameters(self, T1=None, T2=None, M1=None, M2=None,
                            Mt=None, cm_scale=1, silent=False,
                            cm_min=None, cm_max=None, nmt=False,
                            rand_mean=[10, 10, 0.5, 0.5],
                            rand_width=1, save=False):

        if cm_min == None:
            cm_min = self.cm_min
        if cm_max == None:
            cm_max = self.cm_max

        P_ = int(cm_max/(100*cm_scale))
        _P = int(cm_min/(100*cm_scale))

        rand_bounds_ = np.array(rand_mean)*(1 + rand_width)
        _rand_bounds = np.array(rand_mean)*(1 - rand_width)
        _rand_bounds[_rand_bounds < 0] = 0
        rand_bounds_[2:][rand_bounds_[2:] > 1] = 1



        _bounds = [0, 0, 0, 0]
        if Mt == None or nmt:
            bounds_ = [np.inf, np.inf, 1, 1]
        else:
            rand_bounds_[2:][rand_bounds_[2:] > Mt] = Mt
            bounds_ = [np.inf, np.inf, Mt, Mt]

        rand_bounds_ = list(rand_bounds_)
        _rand_bounds = list(_rand_bounds)

        # print(_bounds)
        # print(bounds_)
        # print('rand bound')
        # print(_rand_bounds)
        # print(rand_bounds_)

        variables = [M2, M1, T2, T1]
        fixed = []
        for variable, i in zip(variables, [3, 2, 1, 0]):
            if variable != None:
                bounds_.pop(i)
                _bounds.pop(i)
                _rand_bounds.pop(i)
                rand_bounds_.pop(i)
                fixed.append(i)
        if M1 == None and M2==None and (Mt!=None and not nmt):
            bounds_.pop(-1)
            _bounds.pop(-1)
            _rand_bounds.pop(-1)
            rand_bounds_.pop(-1)
            fixed.append(4)
        variables = [T1, T2, M1, M2]
        T = [T1, T2]
        M = [M1, M2]
        def metr(x):
            j = 0
            for i in range(2):
                if variables[i] == None:
                    T[i] = x[j]
                    j += 1
                else:
                    T[i] = variables[i]

            if M1 == None:
                M[0] = x[j]
                j += 1
            else:
                M[0] = M1
            if M2 != None:
                M[1] = M2
            elif M1 == None and (Mt!=None and not nmt):
                M[1] = (Mt-M[0])/(1-M[0])
            else:
                M[1] = x[j]

            th_E = thld_6(np.array(T, dtype=np.double),
                          np.array(M, dtype=np.double),
                          0, cm_scale, P_)
            return (th_E[_P:P_, _P:P_] - self.WLD.ld[_P:P_, _P:P_]).reshape(-1)


        rng = np.random.default_rng(seed=int(self.WLD.delta*100000%32))
        x0s = rng.uniform(_rand_bounds, rand_bounds_, size=(10, len(rand_bounds_)))
        r = None
        cost = np.inf
        variables_names = ['t1', 't2', 'm1', 'm2', 'mt']
        if fixed != [] and not silent:
            print('fixed variables:' , end=' ')
            for i in fixed:
                print(variables_names[i], end=' ')
        if not silent:
            print('\nfree variables:', end=' ')
            for i in range(len(variables_names)):
                if i not in fixed:
                    print(variables_names[i], end = ' ')
            print()

        res_T = [0, 0]
        res_M = [0, 0]
        for x0 in x0s:
            #print(x0)
            #print(metr(x0))
            res = least_squares(metr, x0, bounds=(_bounds, bounds_))
            if res.cost < cost:
                cost = res.cost
                r = res

        j = 0
        for i in range(2):
            if variables[i] == None:
                res_T[i] = r.x[j]
                j += 1
            else:
                res_T[i] = variables[i]

        if M1 == None:
            res_M[0] = r.x[j]
            j += 1
        else:
            res_M[0] = M1
        if (Mt!=None and not nmt) and (M1 == None and M2 == None):
            res_M[1] = (Mt-res_M[0])/(1-res_M[0])
        elif M2 == None:
            res_M[1] = r.x[j]
        else:
            res_M[1] = M2
        print('estimated:',r.x)
        print('cost function value:', cost)

        if save:
            th_E = thld_6(np.array(res_T, dtype=np.double),
                          np.array(res_M, dtype=np.double),
                          0, cm_scale, P_)
            np.savetxt(self.WLD.ld_filename+'.fld', th_E[_P:P_, _P:P_])
            np.savetxt(self.WLD.ld_filename+'.tld', self.WLD.ld[_P:P_, _P:P_])

        return res_T[0], res_T[1], res_M[0], res_M[1]
