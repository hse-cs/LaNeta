from LaNeta.fastc import binning_b, binning_c, FFT, thld_6, binning_deltas, est_deltas, \
noFFT, binning_b_int, estimate_affine_term, estimate_coef_slow
from LaNeta.fastc import read_positions, read_genotype
import numpy as np
from datetime import datetime
from scipy.optimize import minimize
from scipy.optimize import least_squares

from cyvcf2 import VCF
import LaNeta.gen

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

def read_msprime(ts_list, chr_i, c_max = None, pos_read=True):
    data = Chromosome()
    ts = ts_list[chr_i]
    data.n0 = len(ts.samples(0))
    data.n1 = len(ts.samples(1))
    data.n2 = len(ts.samples(2))

    pos = np.empty(ts.sites().length, dtype=np.double)
    H = np.empty((ts.sites().length, data.n0), dtype=np.double)     # Admixed
    F = np.empty((ts.sites().length, data.n1), dtype=np.double)     # 1st anc
    G = np.empty((ts.sites().length, data.n2), dtype=np.double)    # 2nd anc
    #print(H.shape, F.shape, G.shape)
    c = 0
    if pos_read:
        if data.n1 == 0:
            H = np.empty((ts.sites().length, data.n0 - data.n0//2), dtype=np.double)     # 1st anc
            F = np.empty((ts.sites().length, data.n0//2), dtype=np.double)    # 2nd anc
            data.missing_ref = 1
        if data.n2 == 0:
            H = np.empty((ts.sites().length, data.n0 - data.n0//2), dtype=np.double)     # 1st anc
            G = np.empty((ts.sites().length, data.n0//2), dtype=np.double)    # 2nd anc
            data.missing_ref = 2
        for variant in ts.variants():
            pos[c] = variant.site.position / ts.sequence_length
            if data.n1 == 0:
                F[c] = variant.genotypes[: data.n0 // 2]
                H[c] = variant.genotypes[data.n0 // 2: data.n0]
                G[c] = variant.genotypes[data.n0:]
            elif data.n2 == 0:
                G[c] = variant.genotypes[: data.n0 // 2]
                H[c] = variant.genotypes[data.n0 // 2: data.n0]
                F[c] = variant.genotypes[data.n0:]
            else:
                G[c] = variant.genotypes[data.n0 + data.n1 : data.n0 + data.n1 + data.n2]
                F[c] = variant.genotypes[data.n0 : data.n0 + data.n1]
                H[c] = variant.genotypes[: data.n0]
            c+=1
    else:
        for variant in ts.variants():
            H[c] = variant.genotypes[: data.n0]
            c+=1
            if c == c_max:
                H = H[:c]
                break
    data.G = G
    data.H = H
    data.F = F
    data.n0 = H.shape[1]
    data.n1 = F.shape[1]
    data.n2 = G.shape[1]
    data.pos = pos[:c]
    data.c_x = c
    #print('pa',H.shape)
    return data

def read_pop(samples, popfile, populations):
    H_pop, F_pop, G_pop = populations
    H_samples = []
    F_samples = []
    G_samples = []

    f = open(popfile)
    for line in f:
        sample, pop = line.split()
        sample = sample + '_' + sample #AAAAAAAAAA ANASTASIA IMPROVE IT PLZ!
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

def read_real(files, chr_i, c_max = None, pos_read=True):
    #open vcf
    data = Chromosome()
    vcffile = files[0]
    popfile = files[1]
    mapfile = files[2]
    populations = files[3], files[4], files[5]

    vcf = VCF(vcffile)
    samples = vcf.samples
    H_samples, F_samples, G_samples = read_pop(samples, popfile, populations)
    data.n1 = len(F_samples)
    data.n2 = len(G_samples)
    data.n0 = len(H_samples)
    #read vcf and txt
    F, G, H, pos = None, None, None, None
    if c_max==None:
        c_max = vcf.seqlens[chr_i]
    if pos_read:
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
        pos, c = read_positions(mapfile, c_max, chr_i+1) #chr_i+1 because .txt contains names of chomosomes
        pos = pos[:c]
        H, F, G = read_genotype(vcf(str(chr_i+1)), c, H_samples, F_samples, G_samples)
    else:
        c = c_max
        H, _, _ = read_genotype(vcf(str(chr_i+1)), c, H_samples, F_samples, G_samples)

    data.G = G
    data.H = H
    data.F = F
    #pos = np.linspace(0, 2.86279, c)
    data.pos = pos
    data.c_x = c
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


class Algorithm:
    def __init__(self, exp):
        self.exp = exp

    def estimate_coef_naive(self):
        exp = self.exp
        eps = exp.eps
        self.numerator, self.denumerator = estimate_coef_slow(exp.data.H,   exp.data.n0,
                                              exp.data.F,   exp.data.n1,
                                              exp.data.G,   exp.data.n2,
                                              exp.data.pos, exp.data.c_x,
                                              exp.d_unit,
                                              eps)

        self.numerator, self.denumerator = np.array(self.numerator), np.array(self.denumerator)
        exp.denumerator = self.denumerator
        exp.numerator = self.numerator
        return self.get_coef()


    def estimate_affine_term(self):
        exp = self.exp
        c_chroms = np.zeros((len(exp.chroms), exp.max_c_af, exp.chroms[0].n0), dtype=np.double)
        for i in range(len(exp.chroms)):
            c_chroms[i, :, :] = exp.chroms[i].H

        at = estimate_affine_term(c_chroms, len(exp.chroms), exp.max_c_af, exp.chroms[0].n0)

        return at

    def estimate_bin_b_int(self):
        exp = self.exp
        b = np.array(binning_b_int(exp.data.H,   exp.data.n0,
                                   exp.data.F,   exp.data.n1,
                                   exp.data.G,   exp.data.n2,
                                   exp.c,
                                   exp.data.pos, exp.data.c_x,
                                   exp.d_unit, exp.r,
                                   m=exp.m, missing_ref=exp.data.missing_ref) )
        # print(b.shape)
        # #KILL VARIANCE
        # b[:, self.c > 1000] = 0
        # b[:, self.c < 30] = 0
        exp.b_int.append(b)
        return b

    def estimate_bin_b(self):
        exp = self.exp
        b = np.array(binning_b(exp.data.H,   exp.data.n0,
                               exp.data.F,   exp.data.n1,
                               exp.data.G,   exp.data.n2,
                               exp.data.pos, exp.data.c_x,
                               exp.d_unit) )
        #print(b.shape)
        exp.b.append(b)
        return b

    def estimate_bin_c(self):
        exp = self.exp
        c = np.array(binning_c(exp.data.pos, exp.data.c_x,
                               exp.d_unit, exp.r) )
        exp.c = c
        #print(c)
        i = -1
        while c[i] == 0:
            i-=1
        self.border = len(c) + i + 1
        #print(self.border)
        return c

    def estimate_deltas(self):
        exp = self.exp
        delta = est_deltas(  exp.data.H, exp.data.n0,
                             exp.data.F,   exp.data.n1,
                             exp.data.G,   exp.data.n2,
                             exp.data.pos, exp.data.c_x,
                             m=exp.m, missing_ref=exp.data.missing_ref)
        return delta

    def estimate_deltas_fft(self): #dont use
        print('estimating weights.. ', end='')
        exp = self.exp
        P = exp.P
        def direct():
            deltas_bin = np.array(binning_deltas(exp.data.F,   exp.data.n1,
                                                 exp.data.G,   exp.data.n2,
                                                 exp.data.pos, exp.data.c_x,
                                                 exp.d_unit) )
            deltas_numer = FFT(deltas_bin)
            if self.denumerator == np.empty(0):
                self.denumerator = self.estimate_denumerator
            deltas = np.zeros([P, P])
            for d in range(P):
                for ds in range(P-d):
                    if self.denumerator[d][ds] == 0:
                        deltas[d][ds] += np.nan
                    else:
                        deltas[d][ds] = deltas_numer[d][ds]/self.denumerator[d][ds]
            return deltas
        deltas = direct()
        if exp.unb:
            exp.data.swap_ref()
            deltas += direct()
            deltas /= 2
        print('done!')
        return deltas

    def estimate_numerator(self):
        exp = self.exp
        b = self.estimate_bin_b_int() #no epsilon during binning
        #b = b / (exp.data.n0*exp.data.n1*exp.data.n2) #back to freq
        #print(b[:,:self.border].shape)
        self.numerator = FFT(b[:,:self.border])[1:self.border, 1:self.border]
        #self.numerator = np.round(self.numerator) / ( (exp.data.n0*exp.data.n1*exp.data.n2)**3 )
        #exp.numerator = self.numerator
        return self.numerator

    def estimate_numerator_slow(self):
        exp = self.exp
        b = self.estimate_bin_b_int()
        #print(b[:,:self.border].shape)
        self.numerator = self.numerator / ( (exp.data.n0*exp.data.n1*exp.data.n2)**3 )
        self.numerator = noFFT(b[:,:self.border], b.shape[0], b.shape[1])
        return self.numerator

    def estimate_denumerator_slow(self):
        exp = self.exp
        c = self.estimate_bin_c()
        c = c[:self.border].reshape(1, -1)
        c = np.array(c, dtype=np.double)
        self.numerator = noFFT(c, c.shape[0], c.shape[1])
        return self.numerator

    def estimate_denumerator(self):
        exp = self.exp
        self.c = self.estimate_bin_c()
        c = self.c
        #print(c.shape)
        self.denumerator = FFT(c[:self.border])[1:self.border, 1:self.border] #zeros in the end of array
        self.denumerator = np.round(self.denumerator) #kill machine epsilon (everewhere int)
        return self.denumerator

    def get_coef(self):
        exp = self.exp
        n = exp.data.n0
        #print(n)
        coef = np.divide(exp.numerator, exp.denumerator, out=np.zeros_like(exp.numerator), where=exp.denumerator!=0)*(n/((n-1)*(n-2))) #default out is zeros, where divisor has zero, there is no devision
        return coef

    def estimate_weighted_ld_coef(self):
        print('estimating coef.. ', end='')
        exp = self.exp
        def direct():
            self.denumerator = self.estimate_denumerator()
            self.numerator = self.estimate_numerator()
            coef = self.get_coef()
            return coef
        coef = direct()
        if exp.unb:
            exp.data.swap_ref()
            coef += direct()
            coef /= 2
        print('done!')
        return coef


class ThLd:

    def __init__(self, data_ms=None, mapfile=None, #data_ms - list of ts
                 popfile=None, pop0=None, pop1=None,
                 pop2=None, vcffile=None,
                 T1=10, T2=10, m1=None, m2=None, m=None):
        self.M = [m1, m2]
        self.m = m
        self.m1 = m1
        self.m2 = m2
        if m==None:
            self.m = -1.0
        self.T = [T1, T2]
        self.at = 0
        self.data_ms = data_ms
        self.vcffile = vcffile
        self.mapfile = mapfile
        self.pop0, self.pop1, self.pop2 = pop0, pop1, pop2
        self.popfile = popfile
        self.numerator = np.empty(0)
        self.denumerator = np.empty(0)
        self.coef = np.empty(0)
        self.delta = 0
        self.b_int = []
        self.b = []

    def getCoef(self, cm_min=0.5, cm_max=20):
        P_ = int(cm_max/(100*self.d_unit))
        _P = int(cm_min/(100*self.d_unit))
        return self.coef[_P:P_, _P:P_]

    def getLAcov(self, cm_min=0.5, cm_max=20):
        return self.getCoef(cm_min=cm_min, cm_max=cm_max)/self.delta**3

    def start_slow(self):
            g_num = 0
            for chr_i in self.chroms:
                alg = Algorithm(self)
                self.data = chr_i
                self.denumerator += alg.estimate_denumerator_slow()
                self.numerator += alg.estimate_numerator_slow()
                self.deltas += alg.estimate_deltas()
                g_num+=1
                print('{}/{}'.format(g_num, len(self.chroms)))
            alg = Algorithm(self)
            alg.numerator = self.numerator
            alg.denumerator = self.denumerator
            self.coef = alg.get_coef()/2
            self.deltas /= g_num

    def start_slow_real(self, d_unit=0.001, eps=0.000005):
        self.eps = eps
        self.d_unit = d_unit
        if self.file!=None and self.mapfile!=None:

            vcf = VCF(self.file)
            mapm = open(self.mapfile, 'r')
            g_num = 0

            print('sequence: '+vcf.seqnames[0])
            self.chroms = []
            self.chroms = [gen.realdata(vcf, mapm, 0, self.user_n)]
            self.chr_n = len(self.chroms)
            self.data = self.chroms[0]
            g_num+=1

            alg = Algorithm(self)
            self.coef = alg.estimate_coef_naive()
            self.delta += alg.estimate_deltas()
            print('{}/{}'.format(g_num, len(vcf.seqnames)))

            vcf.close()
            mapm.close()


    def start_affine(self, max_c_af=10000):
        self.max_c_af = max_c_af
        g_num = 0
        self.chroms = []

        if self.data_ms==None:
            vcf = VCF(self.vcffile)
            chr_n = len(vcf.seqlens) # 20
            vcf.close()
            read = read_real
            raw = [self.vcffile, self.popfile, self.mapfile, self.pop0, self.pop1, self.pop2]
        else:
            chr_n =len(self.data_ms)
            read = read_msprime
            raw = self.data_ms

        for i in range(chr_n):
            data = read(raw, i, pos_read=False, c_max = max_c_af)
            self.chroms.append(data)
            g_num+=1
            if g_num == 22:
                break

        alg = Algorithm(self)
        self.at = alg.estimate_affine_term()

    def estimate_time(self, m1=None, m2=None, mt=None, jk=False, du=0.001, r=0.0005, cm_min=1, cm_max=30, af=False):
        if af:
            print('estimating affine term...')
            self.start_affine()
        print('estimating weighted LD...')
        self.start(du=du, r=r)
        if self.m1!=None and self.m2!=None:
            T1, T2 = self.estimate_ls(cm_min=cm_min, cm_max=cm_max)
        else:
            T1, T2 = self.estimate_ls_one_prop(cm_min=cm_min, cm_max=cm_max)
        T1_do = np.zeros((self.chr_n))
        T2_do = np.zeros((self.chr_n))
        T1_ps = np.zeros((self.chr_n))
        T2_ps = np.zeros((self.chr_n))
        _T1, T1_ = np.nan, np.nan
        _T2, T2_ = np.nan, np.nan
        if jk:
            print('running jackknife...')
            for chr_i in range(self.chr_n):
                alg = Algorithm(self)
                self.numerator = np.delete(self.numerators, chr_i, axis=0).sum(axis=0)
                self.denumerator = np.delete(self.denumerators, chr_i, axis=0).sum(axis=0)
                self.numerator /= self.chr_n - 1
                self.denumerator /= self.chr_n - 1
                self.coef = alg.get_coef() - self.at
                self.delta = np.delete(self.deltas, chr_i, axis=0).sum()
                self.delta /= self.chr_n - 1
                if self.m1!=None and self.m2!=None:
                    T1_do[chr_i], T2_do[chr_i] = self.estimate_ls(cm_min=cm_min, cm_max=cm_max)
                else:
                    T1_do[chr_i], T2_do[chr_i] = self.estimate_ls_one_prop(cm_min=cm_min, cm_max=cm_max)
                #print('jk estimation without {} chrom'.format(chr_i+1) , T1_do[chr_i], T2_do[chr_i])
            #print(self.chr_n)
            #print(T1_do)
            #print(T2_do)
            T1_ps = np.abs(self.chr_n * T1 - (self.chr_n-1) * T1_do)
            T2_ps = np.abs(self.chr_n * T2 - (self.chr_n-1) * T2_do)
            #print(T1_ps)
            #print(T2_ps)
            _T1 = T1_ps.mean() - 1.96*np.sqrt((1/self.chr_n)*T1_ps.var())
            T1_ = T1_ps.mean() + 1.96*np.sqrt((1/self.chr_n)*T1_ps.var())
            _T2 = T2_ps.mean() - 1.96*np.sqrt((1/self.chr_n)*T2_ps.var())
            T2_ = T2_ps.mean() + 1.96*np.sqrt((1/self.chr_n)*T2_ps.var())
        print(T1, T2, [_T1, T1_] , [_T2, T2_])
        return T1, T2, [_T1, T1_] , [_T2, T2_]

    def start(self, du=0.001, r=0.0005):

        self.numerator = np.empty(0)
        self.denumerator = np.empty(0)
        self.coef = np.empty(0)
        self.delta = 0
        #data init
        self.d_unit = du
        self.r = r
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
        else:
            chr_n =len(self.data_ms)
            if len(self.data_ms[0].samples(1)) == 0 or len(self.data_ms[0].samples(2)) == 0:
                one_ref = True
            read = read_msprime
            raw = self.data_ms

        #print('CHR_N:', chr_n)
        self.chr_n = chr_n
        g_num = 0 #chr loop

        self.numerators = []
        self.denumerators = []
        self.deltas = []

        for i in range(chr_n):
            delta_i = 0
            self.data = read(raw, i)
            alg = Algorithm(self)
            t0 = datetime.now()
            denumerator_i = alg.estimate_denumerator()
            numerator_i = alg.estimate_numerator()
            delta_i += alg.estimate_deltas()
            if one_ref:
                self.data.swap_ref()
                numerator_i = smart_add(numerator_i, alg.estimate_numerator())
                numerator_i /= 2
            g_num+=1
            #print('delta_i:',delta_i)
            print('{}/{}'.format(g_num, chr_n), datetime.now()-t0)
            self.numerator = smart_add(self.numerator, numerator_i)
            self.denumerator = smart_add(self.denumerator, denumerator_i)
            self.delta += delta_i

            self.numerators.append(numerator_i)
            self.denumerators.append(denumerator_i)
            self.deltas.append(delta_i)

        self.numerators = mkarr(self.numerators)
        self.denumerators = mkarr(self.denumerators)
        self.deltas = np.array(self.deltas)
        #final estimation
        alg = Algorithm(self)
        #print(g_num)
        self.numerator /= g_num
        self.denumerator /= g_num
        self.coef = alg.get_coef() - self.at
        self.delta /= g_num

    def estimate_ls(self, cm_min=1, cm_max=20):
        #print('estimating parameteres ...')

        P_ = int(cm_max/(100*self.d_unit))
        _P = int(cm_min/(100*self.d_unit))

        def metr(x):
            #print(x)
            M = [self.M[0], self.M[1]] #TMNP
            th_E = thld_6(np.array(x, dtype=np.double), np.array(M, dtype=np.double), 0, self.d_unit, P_)
            return (th_E[_P:P_, _P:P_] - self.coef[_P:P_, _P:P_]/self.delta**3).reshape(-1)

        x0 = np.array([10, 10]) # T1, T2
        res = least_squares(metr, x0, bounds=((0, 0), (np.inf, np.inf)))
        #print(res)
        self.T = res.x
        #print(self.T)
        return res.x

    def estimate_ls_one_prop(self, cm_min=1, cm_max=20):
        #print('estimating parameteres ...')

        P_ = int(cm_max/(100*self.d_unit))
        _P = int(cm_min/(100*self.d_unit))

        def metr(x):
            #print(x)
            Mt = self.m #total adm prop
            th_E = thld_6(np.array(x, dtype=np.double), np.array([x[2], (Mt-x[2])/(1-x[2])], dtype=np.double), 0, self.d_unit, P_)
            return (th_E[_P:P_, _P:P_] - self.coef[_P:P_, _P:P_]/self.delta**3).reshape(-1)

        x0 = np.array([10, 10, 0.5]) # T1, T2
        res = least_squares(metr, x0, bounds=((0, 0, 0), (np.inf, np.inf, 1)))
        #print(res)
        self.T = res.x[:2]
        self.M = [res.x[2], (self.m - res.x[2])/(1 - res.x[2])]
        #print(self.T)
        return res.x[:2]

    def estimate(self):
        print('estimating parameteres ...')

        def metr(var, D, a_coef):
            print(var)
            d=5
            P = len(a_coef)
            T = [var[0], var[1]]
            M = [self.M[0], self.M[1]] #TMNP
            th_E = thld_6(np.array(T, dtype=np.double), np.array(M, dtype=np.double), 0, P, self.d_unit)
            return ((th_E[:P//d,:P//d] - a_coef[:P//d,:P//d]/(D**3))**2).sum()

        var = np.array([10, 10]) # T1, T2
        res = minimize(metr, var, (self.deltas, self.coef), method='BFGS')
        print(res)
        self.T = res.x
        print(self.T)
        return res.x


    def estimate_prop(self): #needs update
        print('estimating parameteres ...')

        def metr(var, D, a_coef):
            print(var)
            d=5
            P = len(a_coef)
            T = [self.T[0], self.T[1]]
            M = [var[0], var[1]] #TMNP
            th_E = thld_6(np.array(T, dtype=np.double), np.array(M, dtype=np.double), 0, P, self.d_unit)
            return np.sqrt(((th_E[:P//d,:P//d] - a_coef[:P//d,:P//d]/(D[:P//d,:P//d]))**2).sum())

        var = np.array([self.M[0], self.M[1]]) # T1, T2
        res = minimize(metr, var, (self.deltas, self.coef), method='BFGS')
        print(res)
        self.T = res.x
        print(self.M)
        return res.x

    def estimate_2fix(self): #needs update
        print('estimating parameteres with T2 fix...')

        def metr(T1_var, D, a_coef, aft):
            print(T1_var, self.T[1])
            d=5
            P = len(a_coef)
            T = [T1_var, self.T[1]]
            M = [self.M[0], self.M[1]] #TMNP
            th_E = thld_6(np.array(T, dtype=np.double), np.array(M, dtype=np.double), 0, P, self.d_unit)
            return np.sqrt(((th_E[:P//d,:P//d] - a_coef[:P//d,:P//d]/(D[:P//d,:P//d]) + aft)**2).sum())

        x0 = self.T[0] # T1
        res = least_squares(metr, x0, bounds=((0, 0), (np.inf, np.inf)))
        print(res)
        self.T[0] = res.x
        print(self.T)
        return res.x

    def estimate_1fix(self): #needs update
        print('estimating parameteres with T1 fix...')

        def metr(T2_var, D, a_coef):
            print(self.T[0], T2_var)
            d=5
            P = len(a_coef)
            T = [self.T[0], T2_var]
            M = [self.M[0], self.M[1]] #TMNP
            th_E = thld_6(np.array(T, dtype=np.double), np.array(M, dtype=np.double), 0, P, self.d_unit)
            return np.sqrt(((th_E[:P//d,:P//d] - a_coef[:P//d,:P//d]/(D[:P//d,:P//d]))**2).sum())

        T2_var = self.T[1] # T2
        res = minimize(metr, T2_var, (self.deltas, self.coef), method='BFGS')
        print(res)
        self.T[0] = res.x
        print(self.T)
        return res.x
