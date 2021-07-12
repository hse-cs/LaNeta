from fastc import binning_b, binning_c, FFT, thld_6

import numpy as np
from datetime import datetime
from scipy.optimize import minimize

class ThLd:

    def __init__(self, P, chroms):
        self.chroms = chroms
        self.P = P
        self.chr_n = len(chroms)
        self.M = [0.2, 0.2]
        self.T = [10, 10]


    def cppbin(self):
        chr_n = self.chr_n
        for chr_i in range(chr_n):
            print('\nchrom {}/{}'.format(chr_i + 1, chr_n))
            np.save("F.npy", np.array(self.chroms[chr_i].F, dtype=float))
            np.save("G.npy", np.array(self.chroms[chr_i].G, dtype=float))
            np.save("H.npy", np.array(self.chroms[chr_i].H, dtype=float))
            np.save("pos.npy", self.chroms[chr_i].pos)
            print("files created")
            #cfunction
            #...
            #collecting data
            #b, c = np.load("b.npy"), np.load("c.npy")

    def alg(self, deltas=True, coef=True):

        P = self.P

        nominator = np.zeros((2*P - 1, 2*P - 1))    #from article
        denominator = np.zeros((2*P - 1, 2*P - 1))
        delta = 0
        chr_n = self.chr_n
        for chr_i in range(chr_n):
            print('\nchrom {}/{}'.format(chr_i + 1, chr_n))
            if coef:
                b = np.array(binning_b(self.chroms[chr_i].H,   self.chroms[chr_i].n0,
                                       self.chroms[chr_i].F,   self.chroms[chr_i].n1,
                                       self.chroms[chr_i].G,   self.chroms[chr_i].n2,
                                       self.chroms[chr_i].pos, self.chroms[chr_i].c_x,
                                       self.P) )

                c = np.array(binning_c(self.chroms[chr_i].pos, self.chroms[chr_i].c_x,
                                 self.P) )

                nominator_i, denominator_i = FFT(b), FFT(c)
                nominator += nominator_i
                denominator += denominator_i

            if deltas:
                c = 0
                delta_x = 0
                for x in range(len(self.chroms[chr_i].pos)):
                    delta_x += (self.chroms[chr_i].F[x].sum()/self.chroms[chr_i].n1 - self.chroms[chr_i].G[x].sum()/self.chroms[chr_i].n2)**2
                delta_x /= len(self.chroms[chr_i].pos)
                print('deltas on chrom:',chr_i)
                delta += delta_x

        if coef:
            n = self.chroms[chr_i].n0
            self.coef = np.zeros((P,P))
            for d in range(P):
                for ds in range(P-d):
                    if denominator[d][ds] == 0:
                        self.coef[d][ds] += np.nan
                    else:
                        self.coef[d][ds] = ((n/((n-1)*(n-2)))*nominator[d][ds])/denominator[d][ds]
        if deltas:
            delta /= chr_n

            self.deltas = np.zeros((P,P))
            for d in range(P):
                for ds in range(P-d):
                    self.deltas[d][ds] = delta
            self.deltas = self.deltas**3

    def load(self, deltas='', coef='', file=''):
        if deltas != '':
            self.deltas = np.load(deltas)
            print(deltas, '- succesfully loaded!')
        if coef != '':
            self.coef = np.load(coef)
            print(coef, '- succesfully loaded!')
        if file != '':
            self.file = file

    def estimate(self):
        print('estimating parameteres ...')

        def metr(var, D, a_coef):
            print(var)
            d=5
            P = len(a_coef)
            T = [var[0], var[1]]
            M = [self.M[0], self.M[1]] #TMNP
            th_E = thld_6(np.array(T, dtype=np.double), np.array(M, dtype=np.double), 0, P)
            return np.sqrt(((th_E[:P//d,:P//d] - a_coef[:P//d,:P//d]/(D[:P//d,:P//d]))**2).sum())

        var = np.array([self.T[0], self.T[1]]) # T1, T2
        res = minimize(metr, var, (self.deltas, self.coef), method='BFGS')
        print(res)
        self.T = res.x
        print(self.T)
