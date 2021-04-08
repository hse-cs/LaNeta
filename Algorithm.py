import numpy as np
import conv
import sys
from datetime import datetime

P = 300           # THE ONLY P



def binning(chrom, deltas=False):          # (Chromosome obj, bool for deltas)
    # Part 1 - Binning data

    startTime = datetime.now()
    print('Bin: |', end='')

    d_unit = 1/P     # unit of 'd' in (0,1)-length chromosome

    count = 0
    for x in chrom.pos:
        count+=1

    b = np.zeros((chrom.n0, P))
    if deltas == True:
        b = np.zeros((P))

    c = np.zeros((P))

    d_current = 0 # d_current - position in term of d_unit: {0, 1, ... , d_current, ... , P - 1}
    for x in range(count):
        while chrom.pos[x] >= (d_current + 1) * d_unit:
            d_current += 1

        c[d_current] += 1                           # for each variant in d_current bin c += 1

        f_x = chrom.H[x].sum() / chrom.n0                                 # f_x - impirical allile freq.
        sigma_x = chrom.F[x].sum()/chrom.n1 - chrom.G[x].sum()/chrom.n2               # sigma_x - dif. in allile freq. in two src. pop.

        if deltas == True:
            b[d_current] += sigma_x**2
        else:
            for sample_i in range(chrom.n0):
                b[sample_i][d_current] += sigma_x * (chrom.H[x, sample_i] - f_x)      # computing b_i[d_current] for every i

        if (x+1) % (count//50) == 0:
            print('#', end='')
    endTime = datetime.now()
    print('| 100%, done for', endTime - startTime)
    return b, c    # (c, d)


def FFT(bc):
    # Part 2 - FFT

    startTime = datetime.now()
    fft_b = conv.fft_beauty_conv_up(bc[0])
    endTime = datetime.now()
    print(', nominator done for', endTime - startTime)

    startTime = datetime.now()
    fft_c = conv.fft_beauty_conv_up(bc[1])
    endTime = datetime.now()
    print(', denominator done for', endTime - startTime)
    return fft_b, fft_c

def data_alg(data, deltas=False):

    nominator = np.zeros((2*P - 1, 2*P - 1))    #from article
    denominator = np.zeros((2*P - 1, 2*P - 1))

    # chromosomes loop
    for chrom in data:

        print('\nchrom {}/{}'.format(chrom.name + 1, len(data)))

        nominator_i, denominator_i = FFT(binning(chrom))
        nominator += nominator_i
        denominator += denominator_i

    # Part 3 - Produing return
    n = data[0].n0
    a = np.zeros((P,P))
    for d in range(P):
        for ds in range(P-d):
            if denominator[d][ds] == 0:
                a[d][ds] = np.nan
            else:
                if deltas == False:
                    a[d][ds] = ((n/((n-1)*(n-2)))*nominator[d][ds])/denominator[d][ds]
                else:
                    a[d][ds] = nominator[d][ds]/denominator[d][ds]
    return a

def deltas_alg_lob(data): #deltas without d and ds dependence
    delta = 0
    for chrom in data:
        c = 0
        delta_x = 0
        for x in range(len(chrom.pos)):
            delta_x += (chrom.F[x].sum()/chrom.n1 - chrom.G[x].sum()/chrom.n2)**2
        delta_x /= len(chrom.pos)
        delta += delta_x
        print('deltas on chrom:',chrom)
    delta /= len(data)

    deltas_coef = np.zeros((P,P))
    for d in range(P):
        for ds in range(P-d):
            deltas_coef[d][ds] = delta
    return deltas_coef**3



def alg(data, coef=True, deltas=True):
    c = np.empty((0))
    d = np.empty((0))

    if coef == True:
        c = data_alg(data)
    if deltas == True:
        #d = deltas_alg(data)
        d = deltas_alg_lob(data)

    return c, d
