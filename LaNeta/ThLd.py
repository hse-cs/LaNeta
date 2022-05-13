from LaNeta.fastc import binning_b, binning_c, FFT, thld_6, binning_deltas, est_deltas, \
noFFT, binning_b_int, estimate_affine_term, estimate_coef_slow, estimate_affine_term_v2
from LaNeta.fastc import read_vcf_pos, read_vcf_genotype, read_positions
import numpy as np
from datetime import datetime
from scipy.optimize import minimize
from scipy.optimize import least_squares

from cyvcf2 import VCF
import LaNeta.gen


def get_seqlen(vcffile, name):
    vcf = VCF(vcffile)
    c = 0
    for var in vcf(name):
        c += 1
    vcf.close()
    return c

def read_pos(mapfile, c_max, chr_name):
    position, c = read_positions(mapfile, c_max, chr_name)
    pos = np.array(position)
    pos = np.sort(pos)
    print(pos)
    return pos, c

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

def read_msprime(ts_list, chr_i, c_max = None, pos_read=True, gt=False):
    data = Chromosome()
    ts = ts_list[chr_i]
    data.n0 = len(ts.samples(0))# * 10
    data.n1 = len(ts.samples(1))
    data.n2 = len(ts.samples(2))


    pos = np.empty(ts.sites().length, dtype=np.double)
    H = np.empty((ts.sites().length, data.n0), dtype=np.double)     # Admixed
    F = np.empty((ts.sites().length, data.n1), dtype=np.double)     # 1st anc
    G = np.empty((ts.sites().length, data.n2), dtype=np.double)    # 2nd anc
    #print(H.shape, F.shape, G.shape)
    c = ts.sites().length
    if data.n1 == 0:
        H = np.empty((ts.sites().length, data.n0 - data.n0//2), dtype=np.double)     # 1st anc
        F = np.empty((ts.sites().length, data.n0//2), dtype=np.double)    # 2nd anc
        data.missing_ref = 1
    if data.n2 == 0:
        H = np.empty((ts.sites().length, data.n0 - data.n0//2), dtype=np.double)     # 1st anc
        G = np.empty((ts.sites().length, data.n0//2), dtype=np.double)    # 2nd anc
        data.missing_ref = 2
    if pos_read:
        pos = ts.tables.sites.position / ts.sequence_length
        gt_matrix = ts.genotype_matrix()
        if data.n1 == 0:
            F = gt_matrix[:, : data.n0 // 2]
            H = gt_matrix[:, data.n0 // 2: data.n0]
            G = gt_matrix[:, data.n0:]
        elif data.n2 == 0:
            G = gt_matrix[:, : data.n0 // 2]
            H = gt_matrix[:, data.n0 // 2: data.n0]
            F = gt_matrix[:, data.n0:]
        else:
            G = gt_matrix[:, data.n0 + data.n1 : data.n0 + data.n1 + data.n2]
            F = gt_matrix[:, data.n0 : data.n0 + data.n1]
            H = gt_matrix[:, : data.n0]
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
    print(data.n0, data.n1, data.n2, H.dtype)
    data.G = np.array(G[:, :data.n2])
    data.H = np.array(H[:, :data.n0])
    data.F = np.array(F[:, :data.n1])
    data.pos = pos[:c]
    data.c_x = c
    #print('pa',H.shape)

    # if gt:
    #     data.n0*=2
    #     data.n1*=2
    #     data.n2*=2
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
        elif pop == F_pop:
            if sample in samples:
                F_samples.append(samples.index(sample))
        elif pop == G_pop:
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
    if c_max==None:
        c_max = get_seqlen(vcffile ,chr_name)
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
        pos, _ = read_pos(mapfile, c_max, chr_name) #chr_i+1 because .txt contains names of chomosomes
        H, F, G = read_vcf_genotype(vcf(chr_name), c_max, H_samples, F_samples, G_samples)
    else:
        _, c_total = read_pos(mapfile, c_max, chr_name)
        print('max:', c_max, 'total:', c_total)
        rng = np.random.default_rng()
        inds = np.sort(rng.choice(c_total, size=c_max, replace=False))
        H, F, G = read_vcf(vcf(chr_name), c_max, H_samples, F_samples, G_samples, inds=inds)
        c = c_max
    vcf.close()

    data.G = G
    data.H = H
    data.F = F
    #pos = np.linspace(0, 2.86279, c)
    print(np.array(data.H))
    print(np.array(pos))
    data.pos = pos
    data.c_x = c_max
    return data


class Chromosome:
    def __init__(self): # data = [H_gen, F_gen, G_gen, seqlens]
        self.H_vcf, self.F_vcf, self.G_vcf = None, None, None
        self.n0, self.n1, self.n2 = 0, 0, 0
        self.to_swap = None
        self.vcffiles = None
        self.missing_ref = 0

    def swap_ref(self):
        if self.missing_ref == 2:
            self.H, self.G = self.G, self.H
            self.n0, self.n2 = self.n2, self.n0
        elif self.missing_ref == 1:
            self.H, self.F = self.F, self.H
            self.n0, self.n1 = self.n1, self.n0
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
                                              WLD.data.pos, WLD.data.c_x,
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

    def estimate_bin_b_int(self):
        WLD = self.WLD
        print(WLD.data.missing_ref)
        b = np.array(binning_b_int(WLD.data.H,   WLD.data.n0,
                                   WLD.data.F,   WLD.data.n1,
                                   WLD.data.G,   WLD.data.n2,
                                   WLD.freq_filter, WLD.c,
                                   WLD.data.pos, WLD.data.c_x,
                                   WLD.bin_size, WLD.bin_radius,
                                   m=WLD.Mt, missing_ref=WLD.data.missing_ref) )
        #WLD.b_int.append(b)
        return b

    def estimate_bin_b(self):
        WLD = self.WLD
        b = np.array(binning_b(WLD.data.H,   WLD.data.n0,
                               WLD.data.F,   WLD.data.n1,
                               WLD.data.G,   WLD.data.n2,
                               WLD.data.pos, WLD.data.c_x,
                               WLD.bin_size) )
        #print(b.shape)
        WLD.b.append(b)
        return b

    def estimate_bin_c(self):
        WLD = self.WLD
        c = np.array(binning_c(WLD.data.pos, WLD.data.c_x,
                               WLD.bin_size, WLD.bin_radius) )
        WLD.c = c
        #print(c)
        i = -1
        while c[i] == 0:
            i-=1
        self.border = len(c) + i + 1
        # print('pos len', len(WLD.data.pos))
        # print('c_x', WLD.data.c_x)
        # print('border =', self.border)
        return c

    def estimate_deltas(self):
        WLD = self.WLD
        delta = est_deltas(  WLD.data.H, WLD.data.n0,
                             WLD.data.F,   WLD.data.n1,
                             WLD.data.G,   WLD.data.n2,
                             WLD.freq_filter, WLD.data.c_x,
                             m=WLD.Mt, missing_ref=WLD.data.missing_ref)
        return delta

    def estimate_deltas_numer(self): #dont use
        print('estimating weights (FFT).. ')
        WLD = self.WLD
        deltas_bin = np.array(binning_deltas(WLD.data.H,   WLD.data.n0,
                                             WLD.data.F,   WLD.data.n1,
                                             WLD.data.G,   WLD.data.n2,
                                             WLD.data.pos, WLD.data.c_x,
                                             WLD.bin_size,
                                             m=WLD.m, missing_ref=WLD.data.missing_ref) )
        deltas_numer = FFT(deltas_bin[:self.border])[1:self.border, 1:self.border]
        return deltas_numer

    def estimate_numerator(self):
        WLD = self.WLD
        b = self.estimate_bin_b_int() #no epsilon during binning
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
        self.c = self.estimate_bin_c()
        c = self.c
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
                 gt=False):

        self.coef = None
        self.delta = None
        self.ld = None

        #data load init
        if data_ms or vcffile:
            self.one_ref=False
            self.gt = gt
            if data_ms==None:
                vcf = VCF(vcffile)
                self.chr_names = vcf.seqnames # 20
                self.chr_n = len(self.chr_names)
                vcf.close()
                if pop1 == None or pop2 == None:
                        self.one_ref = True
                self.read = read_real
                self.data_input = [vcffile, popfile, mapfile,
                                   pop0, pop1, pop2]
            else:
                self.chr_n = len(data_ms)
                self.chr_names = range(self.chr_n)
                if len(data_ms[0].samples(1)) == 0 or len(data_ms[0].samples(2)) == 0:
                    if self.Mt:
                        self.one_ref = True
                    else:
                        pass
                self.read = read_msprime
                self.data_input = data_ms

    def estimate(self, bin_size=0.001, bin_radius=0.0005, freq_filter=10,
                 at=False, Mt=-1):
        self.Mt = Mt
        self.freq_filter = freq_filter/100
        self.coef = np.empty(0)
        #self.delta_numer = np.empty(0)
        self.bin_size = bin_size
        self.bin_radius = bin_radius
        print('CHR_N:', self.chr_n)
        g_num = 0 #chr loop

        self.numerators = [0]*self.chr_n
        self.denumerators = [0]*self.chr_n
        self.deltas = [0]*self.chr_n
        #self.delta_numers = []

        for chr_name in self.chr_names:
            t0 = datetime.now()
            self.data = self.read(self.data_input, chr_name, gt=self.gt)
            print(f'{datetime.now() - t0} - data read.')
            t0 = datetime.now()
            alg = Weighted_LD_Algorithm(self)
            denumerator_i = alg.estimate_denumerator()
            numerator_i = alg.estimate_numerator()
            #delta_numer_i = alg.estimate_deltas_numer()
            delta_i = alg.estimate_deltas()
            if self.one_ref:
                self.data.swap_ref()
                numerator_i += alg.estimate_numerator()
                numerator_i /= 2

                #delta_numer_i = smart_add(delta_numer_i, alg.estimate_deltas_numer())
                #delta_numer_i /= 2

                delta_i += alg.estimate_deltas()
                delta_i /= 2
            #print('{}/{}'.format(g_num, chr_n), datetime.now()-t0, 'delta_numer_mean:', delta_numer_i.mean())
            #self.delta_numer = smart_add(self.delta_numer, alg.estimate_deltas_numer())
            print(f'{datetime.now() - t0} - estimation.')
            t0 = datetime.now()
            self.numerators[g_num] = numerator_i
            self.denumerators[g_num] = denumerator_i
            #self.delta_numers.append(delta_numer_i)
            self.deltas[g_num] = delta_i
            g_num+=1
            print(f'{datetime.now() - t0} - enter.')
            print('{}/{}'.format(g_num, self.chr_n), 'delta_i:', delta_i)
            #break

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

    def estimate_affine_term(self):
        pass

    def jackknife(self, chr_i):
        alg = Algorithm(self)
        self.numerator = np.delete(self.numerators, chr_i, axis=0).sum(axis=0)
        self.denumerator = np.delete(self.denumerators, chr_i, axis=0).sum(axis=0)
        self.numerator /= self.chr_n - 1
        self.denumerator /= self.chr_n - 1
        self.coef = alg.get_coef() - self.at
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

    def __init__(self, data_ms=None, mapfile=None, #data_ms - list of ts
                 popfile=None, pop0=None, pop1=None,
                 pop2=None, vcffile=None, gt=False):
        self.Ms = []
        self.WLD = Weighted_LD(data_ms=data_ms, mapfile=mapfile, vcffile=vcffile,
                               popfile=popfile, pop0=pop0, pop1=pop1, pop2=pop2,
                               gt=gt)

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
            vcf_len = get_seqlen(self.vcffile, chr_name)
            print('vcf len:',vcf_len)
            pos, txt_len = read_pos(self.mapfile, vcf_len*2, chr_name)
            print('txt len:',txt_len)
            if vcf_len != txt_len:
                print('number of variants in .txt and .vcf is different!')



        #check pop

        #check calculations

    def getwld(self, cm_min=0.5, cm_max=20):
        return self.WLD.getwld(cm_min=cm_min, cm_max=cm_max)

    def getld(self, cm_min=0.5, cm_max=20):
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
                               cm_min=1, cm_max=30, at=False, Mt=-1,
                               freq_filter=10,
                               fixed_time=0):
        print('estimating weighted LD...')
        self.WLD.estimate(bin_size=bin_size, bin_radius=bin_radius,
                          freq_filter=freq_filter, Mt=Mt, at=at)
        parameters = self.estimate_parameters(T1=T1, T2=T2, M1=M1, M2=M2,
                                              Mt=Mt,
                                              cm_scale=bin_size,
                                              cm_min=cm_min, cm_max=cm_max)
        if jk:
            parameters_do = np.zeros((self.WLD.chr_n, parameters.shape[0]))
            _parameters = np.zeros(parameters.shape[0])
            parameters_ = np.zeros(parameters.shape[0])
            print('running jackknife...')
            for chr_i in range(self.WLD.chr_n):
                self.WLD.jackknife(chr_i)
                parameters_do[chr_i] = self.estimate_parameters(T1=T1, T2=T2, M1=M1, M2=M2,
                                                                Mt=self.WLD.Mt,
                                                                cm_scale=self.WLF.bin_size,
                                                                cm_min=0.5, cm_max=20)
            parameters_ps = np.abs(self.WLD.chr_n * parameters - (self.WLD.chr_n-1) * parameters_do)
            # _parameters = parameters_do.mean(axis=0) - 1.96*np.sqrt((1/self.WLD.chr_n)*parameters_do.var(axis=0))
            # parameters_ = parameters_do.mean(axis=0) + 1.96*np.sqrt((1/self.WLD.chr_n)*parameters_do.var(axis=0))
            _parameters = parameters_ps.mean(axis=0) - 1.96*np.sqrt((1/self.WLD.chr_n)*parameters_ps.var(axis=0))
            parameters_ = parameters_ps.mean(axis=0) + 1.96*np.sqrt((1/self.WLD.chr_n)*parameters_ps.var(axis=0))
            return parameters, _parameters, parameters_
        return parameters

    def estimate_parameters(self, T1=None, T2=None, M1=None, M2=None,
                            Mt=None, cm_scale=1,
                            cm_min=0, cm_max=50):

        P_ = int(cm_max/(100*cm_scale))
        _P = int(cm_min/(100*cm_scale))

        rand_bounds_ = [10, 10, 1, 1]
        _rand_bounds = [0, 0, 0, 0]
        _bounds = [0, 0, 0, 0]
        bounds_ = [np.inf, np.inf, 1, 1]
        variables = [M2, M1, T2, T1]
        fixed = []
        for variable, i in zip(variables, [3, 2, 1, 0]):
            if variable != None:
                bounds_.pop(i)
                _bounds.pop(i)
                _rand_bounds.pop(i)
                rand_bounds_.pop(i)
                fixed.append(i)
        if M1 == None and M2==None and Mt!=-1:
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
            if Mt != -1 and (M1 == None and M2 == None):
                M[1] = (Mt-M[0])/(1-M[0])
            elif M2 == None:
                M[1] = x[j]
            else:
                M[1] = M2
            th_E = thld_6(np.array(T, dtype=np.double),
                          np.array(M, dtype=np.double),
                          0, cm_scale, P_)
            return (th_E[_P:P_, _P:P_] - self.WLD.ld[_P:P_, _P:P_]).reshape(-1)


        rng = np.random.default_rng(seed=int(self.WLD.delta*100000%32))
        x0s = rng.uniform(_rand_bounds, rand_bounds_, size=(10, len(rand_bounds_)))
        r = None
        cost = 100.0
        variables_names = ['T1', 'T2', 'M1', 'M2', 'Mt']
        if fixed != []:
            print('fixed variables:' , end=' ')
            for i in fixed:
                print(variables_names[i], end=' ')
            print()

        for x0 in x0s:
            print(x0)
            print(metr(x0))
            res = least_squares(metr, x0, bounds=(_bounds, bounds_))
            print(res.x, res.cost)
            if res.cost < cost:
                cost = res.cost
                r = res
        print('estimated (rand.):',r.x, cost)
        return T[0], T[1], M[0], M[1]
