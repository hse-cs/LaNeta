from fastc import binning_b, binning_c, FFT, thld_6, binning_deltas, est_deltas, \
noFFT, binning_b_int, estimate_affine_term, estimate_coef_slow
from fastc import read_positions, read_genotype
import numpy as np
from datetime import datetime
from scipy.optimize import minimize
from scipy.optimize import least_squares

from cyvcf2 import VCF
import gen

def smart_add(a, b):
    if a.shape[0]>b.shape[0]:
        a[:b.shape[0], :b.shape[0]] += b
        return a
    else:
        b[:a.shape[0], :a.shape[0]] += a
        return b

class Chromosome:
    def __init__(self): # data = [H_gen, F_gen, G_gen, seqlens]
        self.H_vcf, self.F_vcf, self.G_vcf = None, None, None
        self.n0, self.n1, self.n2 = 0, 0, 0
        self.vcffiles = None

    def read_real(self, files, chr_i, c_max = None, pos_read=True):
        #open vcf
        vcffiles = files[0]
        mapfile = files[1]
        H_vcf = VCF(vcffiles[0])
        if vcffiles[1] != None:
            F_vcf = VCF(vcffiles[1])
            self.n1 = len(F_vcf.samples)
        elif vcffiles[2] != None:
            G_vcf = VCF(vcffiles[2])
            self.n2 = len(G_vcf.samples)

        self.n0  = len(H_vcf.samples)
        #read vcf and txt

        F, G, H, pos = None, None, None, None
        if c_max==None:
            c_max = H_vcf.seqlens[chr_i]
        if pos_read:
            pos, c = read_positions(mapfile, c_max, chr_i+1) #chr_i+1 because .txt contains names of chomosomes
            pos = pos[:c]
            if vcffiles[1] == None:
                H = read_genotype(H_vcf(str(chr_i+1)), c, n_min=0, n_max=self.n0//2)
                F = read_genotype(H_vcf(str(chr_i+1)), c, n_min=self.n0//2, n_max=self.n0)
                G = read_genotype(G_vcf(str(chr_i+1)), c, n_max=self.n2)
                self.n1 = self.n0 - self.n0//2
                self.n0 = self.n0//2
            elif vcffiles[2] == None:
                H = read_genotype(H_vcf(str(chr_i+1)), c, n_min=0, n_max=self.n0//2)
                G = read_genotype(H_vcf(str(chr_i+1)), c, n_min=self.n0//2, n_max=self.n0)
                F = read_genotype(F_vcf(str(chr_i+1)), c, n_max=self.n1)
                self.n2 = self.n0 - self.n0//2
                self.n0 = self.n0//2
            else:
                H = read_genotype(H_vcf(str(chr_i+1)), c, n_max=self.n0)
                F = read_genotype(F_vcf(str(chr_i+1)), c, n_max=self.n1)
                G = read_genotype(G_vcf(str(chr_i+1)), c, n_max=self.n2)
        else:
            c = c_max
            H = read_genotype(H_vcf(str(chr_i+1)), c, n_max=self.n0)

        self.G = G
        self.H = H
        self.F = F
        #pos = np.linspace(0, 2.86279, c)
        self.pos = pos
        self.c_x = c

    def read_msprime(self, ts_list, chr_i, c_max = None, pos_read=True):
        ts = ts_list[chr_i]
        self.n0 = len(ts.samples(0))
        self.n1 = len(ts.samples(1))
        self.n2 = len(ts.samples(2))

        pos = np.empty(ts.sites().length, dtype=np.double)
        H = np.empty((ts.sites().length, self.n0), dtype=np.byte)     # Admixed
        F = np.empty((ts.sites().length, self.n1), dtype=np.byte)     # 1st anc
        G = np.empty((ts.sites().length, self.n2), dtype=np.byte)    # 2nd anc
        if self.n1 == 0:
            H = np.empty((ts.sites().length, self.n0 - self.n0//2), dtype=np.byte)     # 1st anc
            F = np.empty((ts.sites().length, self.n0//2), dtype=np.byte)    # 2nd anc
        if self.n2 == 0:
            H = np.empty((ts.sites().length, self.n0 - self.n0//2), dtype=np.byte)     # 1st anc
            G = np.empty((ts.sites().length, self.n0//2), dtype=np.byte)    # 2nd anc
        #print(H.shape, F.shape, G.shape)
        c = 0
        if pos_read:
            for variant in ts.variants():
                pos[c] = variant.site.position / ts.sequence_length
                if self.n1 == 0:
                    F[c] = variant.genotypes[: self.n0 // 2]
                    H[c] = variant.genotypes[self.n0 // 2: self.n0]
                    G[c] = variant.genotypes[self.n0:]
                    self.n1 = self.n0 - self.n0//2
                    self.n0 = self.n0//2
                elif self.n2 == 0:
                    G[c] = variant.genotypes[: self.n0 // 2]
                    H[c] = variant.genotypes[self.n0 // 2: self.n0]
                    F[c] = variant.genotypes[self.n0:]
                    self.n2 = self.n0 - self.n0//2
                    self.n0 = self.n0//2
                else:
                    G[c] = variant.genotypes[self.n0 + self.n1 : self.n0 + self.n1 + self.n2]
                    F[c] = variant.genotypes[self.n0: self.n0 + self.n1]
                    H[c] = variant.genotypes[: self.n0]
                c+=1
        else:
            for variant in ts.variants():
                H[c] = variant.genotypes[: self.n0]
                c+=1
                if c == c_max:
                    H = H[:c]
                    break
        self.G = G
        self.H = H
        self.F = F
        self.pos = pos[:c]
        self.c_x = c

    def swap_ref(self):
        self.H, self.G = self.G, self.H
        self.n0, self.n2 = self.n2, self.n0
        print('ref and adm swaped!')


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
        c_chroms = np.zeros((len(exp.chroms), exp.max_c_af, exp.chroms[0].n0), dtype=np.byte)
        for i in range(len(exp.chroms)):
            c_chroms[i, :, :] = exp.chroms[i].H

        at = estimate_affine_term(c_chroms, len(exp.chroms), exp.max_c_af, exp.chroms[0].n0)

        return at

    def estimate_bin_b_int(self):
        exp = self.exp
        b = np.array(binning_b_int(exp.data.H,   exp.data.n0,
                                   exp.data.F,   exp.data.n1,
                                   exp.data.G,   exp.data.n2,
                                   exp.data.pos, exp.data.c_x,
                                   exp.d_unit) )
        exp.b_int.append(b)
        return b

    def estimate_bin_b(self):
        exp = self.exp
        b = np.array(binning_b(exp.data.H,   exp.data.n0,
                               exp.data.F,   exp.data.n1,
                               exp.data.G,   exp.data.n2,
                               exp.data.pos, exp.data.c_x,
                               exp.d_unit) )
        print(b.shape)
        exp.b.append(b)
        return b

    def estimate_bin_c(self):
        exp = self.exp
        c = np.array(binning_c(exp.data.pos, exp.data.c_x,
                               exp.d_unit) )
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
        delta = est_deltas(exp.data.F,   exp.data.n1,
                             exp.data.G,   exp.data.n2,
                             exp.data.pos, exp.data.c_x)
        delta = delta**3
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
        b = b / (exp.data.n0*exp.data.n1*exp.data.n2) #back to freq
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
        c = self.estimate_bin_c()
        self.denumerator = FFT(c[:self.border])[1:self.border, 1:self.border] #zeros in the end of array
        self.denumerator = np.round(self.denumerator) #kill machine epsilon (everewhere int)
        return self.denumerator

    def get_coef(self):
        exp = self.exp
        n = exp.data.n0
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
                 vcf0=None, vcf1=None, vcf2=None,
                 T1=10, T2=10, m1=0.2, m2=0.2):
        self.M = [m1, m2]
        self.T = [T1, T2]
        self.at = 0
        self.data_ms = data_ms
        self.vcffiles = [vcf0, vcf1, vcf2]
        self.mapfile = mapfile
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
        if self.vcffiles[0] != None:
            vcf = VCF(self.vcffiles[0])
            for i in range(len(vcf.seqnames)):
                data = Chromosome()
                data.read_real([self.vcffiles, self.mapfile], i, c_max=self.max_c_af, pos_read=False)
                self.chroms.append(data)
                self.chr_n = len(self.chroms)
                g_num+=1
                if g_num == 22:
                    break

            alg = Algorithm(self)
            self.at = alg.estimate_affine_term()

            vcf.close()
        else:
            for i in range(len(self.data_ms)):
                data = Chromosome()
                data.read_msprime(self.data_ms, i, c_max=self.max_c_af, pos_read=False)
                self.chroms.append(data)
                self.chr_n = len(self.chroms)
                g_num+=1
                if g_num == 22:
                    break

            alg = Algorithm(self)
            self.at = alg.estimate_affine_term()

    def start(self, du=0.001):

        self.numerator = np.empty(0)
        self.denumerator = np.empty(0)
        self.coef = np.empty(0)
        self.delta = 0
        #data init
        self.d_unit = du
        self.data = Chromosome()
        one_ref=False
        if self.data_ms==None:
            vcf = VCF(self.vcffiles[0])
            chr_n = len(vcf.seqlens) # 20
            vcf.close()
            if self.vcffiles[1] == None or self.vcffiles[2] == None:
                one_ref = True
            self.data.read = self.data.read_real
            raw = [self.vcffiles, self.mapfile]
        else:
            chr_n =len(self.data_ms)
            if len(self.data_ms[0].samples(1)) == 0 or len(self.data_ms[0].samples(2)) == 0:
                one_ref = True
            self.data.read = self.data.read_msprime
            raw = self.data_ms

        print(chr_n, one_ref)

        g_num = 0 #chr loop
        for i in range(chr_n):
            delta_i = 0
            self.data.read(raw, i)
            alg = Algorithm(self)
            t0 = datetime.now()
            denumerator_i = alg.estimate_denumerator()
            numerator_i = alg.estimate_numerator()
            delta_i += alg.estimate_deltas()
            if one_ref:
                self.data.swap_ref()
                numerator_i = smart_add(numerator_i, alg.estimate_numerator())
                #numerator_i /= 2
                delta_i += alg.estimate_deltas()
                #delta_i /= 2
            g_num+=1
            print('{}/{}'.format(g_num, chr_n), datetime.now()-t0)
            self.numerator = smart_add(self.numerator, numerator_i)
            self.denumerator = smart_add(self.denumerator, denumerator_i)
            self.delta += delta_i

        #final estimation
        alg = Algorithm(self)
        self.numerator /= g_num
        self.denumerator /= g_num
        self.coef = alg.get_coef() - self.at
        self.delta /= g_num

    def estimate_ls(self, cm=20):
        #print('estimating parameteres ...')

        self.P_ = int(cm/(100*self.d_unit))

        def metr(x):
            #print(x)
            M = [self.M[0], self.M[1]] #TMNP
            th_E = thld_6(np.array(x, dtype=np.double), np.array(M, dtype=np.double), 0, self.d_unit, self.P_)
            return (th_E[:self.P_,:self.P_] - self.coef[:self.P_,:self.P_]/self.delta).reshape(-1)

        x0 = np.array([10, 10]) # T1, T2
        res = least_squares(metr, x0, bounds=((0, 0), (np.inf, np.inf)))
        #print(res)
        self.T = res.x
        #print(self.T)
        return res.x

    def estimate_ls_one_prop(self, cm=20):
        #print('estimating parameteres ...')

        self.P_ = int(cm/(100*self.d_unit))

        def metr(x):
            #print(x)
            Mt = self.Mt #total adm prop
            th_E = thld_6(np.array(x, dtype=np.double), np.array([x[2], (Mt-x[2])/(1-x[2])], dtype=np.double), 0, self.d_unit, self.P_)
            return (th_E[:self.P_,:self.P_] - self.coef[:self.P_,:self.P_]/self.delta).reshape(-1)

        x0 = np.array([10, 10, 0.5]) # T1, T2
        res = least_squares(metr, x0, bounds=((0, 0, 0), (np.inf, np.inf, 1)))
        #print(res)
        self.T = res.x[:2]
        self.M = [res.x[2], (self.Mt - res.x[2])/(1 - res.x[2])]
        #print(self.T)
        return res.x

    def estimate(self):
        print('estimating parameteres ...')

        def metr(var, D, a_coef):
            print(var)
            d=5
            P = len(a_coef)
            T = [var[0], var[1]]
            M = [self.M[0], self.M[1]] #TMNP
            th_E = thld_6(np.array(T, dtype=np.double), np.array(M, dtype=np.double), 0, P, self.d_unit)
            return ((th_E[:P//d,:P//d] - a_coef[:P//d,:P//d]/(D[:P//d,:P//d]))**2).sum()

        var = np.array([10, 10]) # T1, T2
        res = minimize(metr, var, (self.deltas, self.coef), method='BFGS')
        print(res)
        self.T = res.x
        print(self.T)
        return res.x


    def estimate_prop(self):
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

    def estimate_2fix(self):
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

    def estimate_1fix(self):
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
