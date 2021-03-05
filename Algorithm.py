import numpy as np
import conv
import sys
from datetime import datetime

P = 300           # THE ONLY P


def handingfunction(chrom):
    b, c = FFT(binning(chrom))
    return b, c


def binning(chrom):          # (Chromosome obj)
    # Part 1 - Binning data
    startTime = datetime.now()
    print('Bin: |', end='')

    d_unit = 1/P     # unit of 'd' in (0,1)-length chromosome

    count = 0
    for x in chrom.pos:
        count+=1

    b = np.zeros((chrom.n0, P))
    c = np.zeros((P))

    d_current = 0 # d_current - position in term of d_unit: {0, 1, ... , d_current, ... , P - 1}
    for x in range(count):
        while chrom.pos[x] >= (d_current + 1) * d_unit:
            d_current += 1

        c[d_current] += 1                           # for each variant in d_current bin c += 1

        f_x = chrom.H[x].sum() / chrom.n0                                 # f_x - impirical allile freq.
        sigma_x = chrom.F[x].sum()/chrom.n1 - chrom.G[x].sum()/chrom.n2               # sigma_x - dif. in allile freq. in two src. pop.

        for sample_i in range(chrom.n0):
            b[sample_i][d_current] += sigma_x * (chrom.H[x, sample_i] - f_x)      # computing b_i[d_current] for every i

        if x % (count//50) == 0:
            print('#', end='')
    endTime = datetime.now()
    print('| 100%, done for', endTime - startTime)
    return b, c    # (c, d)

def FFT(bc):
    startTime = datetime.now()
    # Part 2 - FFT
    fft_b = conv.fft_beauty_conv_up(bc[0])
    fft_c = conv.fft_beauty_conv(bc[1])
    endTime = datetime.now()
    print(', done for', endTime - startTime)
    return fft_b, fft_c

def alg(data):

    nominator = np.zeros((2*P - 1, 2*P - 1))
    denominator = np.zeros((2*P - 1, 2*P - 1))

    # chromosomes loop
    for chrom in data:

        print('\nchrom {}/{}'.format(chrom.name + 1, len(data)))

        nominator_i, denominator_i = handingfunction(chrom)
        nominator += nominator_i
        denominator += denominator_i

    # Part 3 - Produing return
    n = data[0].n0
    w3pld_COEF = np.zeros((P,P))
    for d in range(P):
        for ds in range(P-d):
            if denominator[d][ds] == 0:
                w3pld_COEF[d][ds] = np.nan
            else:
                w3pld_COEF[d][ds] = ((n/((n-1)*(n-2)))*nominator[d][ds])/denominator[d][ds]
    return w3pld_COEF
