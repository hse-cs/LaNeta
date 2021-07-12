# cython: language_level=3

cimport cython
from libc.math cimport exp
import numpy as np
from datetime import datetime
from numpy.fft import fft, ifft2



@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def binning_b(char[:,::1] H, int n0,
              char[:,::1] F, int n1,
              char[:,::1] G, int n2,
              double[:] pos, int c_x,
                             int P):          # (Chromosome obj, bool for deltas)
    # Part 1 - Binning data
    startTime = datetime.now()
    print('Bin: |', end='')

    cdef double d_unit = 1.0/P     # unit of 'd' in (0,1)-length chromosome

    cdef double[:,::1] b = np.zeros((n0, P), dtype=np.double)

    cdef int d_current = 0 # d_current - position in term of d_unit: {0, 1, ... , d_current, ... , P - 1}



    cdef double f_x
    cdef double F_x
    cdef double G_x
    cdef double sigma_x
    cdef int x
    cdef int s_i
    for x in range(c_x):
        while pos[x] >= (d_current + 1) * d_unit:
            d_current += 1

        f_x = 0                                # f_x - impirical allile freq.
        for s_i in range(n0):
            f_x += H[x][s_i]
        f_x = f_x / n0


        F_x = 0                                # f_x - impirical allile freq.
        for s_i in range(n1):
            F_x += F[x][s_i]
        F_x = F_x / n1

        G_x = 0                                # f_x - impirical allile freq.
        for s_i in range(n2):
            G_x += G[x][s_i]
        G_x = G_x / n2

        sigma_x = (F_x - G_x) # sigma_x - dif. in allile freq. in two src. pop.

        #if deltas == True:
        #    self.b[d_current] += sigma_x**2
        #else:
        for s_i in range(n0):
            b[s_i, d_current] += sigma_x * (H[x, s_i] - f_x)      # computing b_i[d_current] for every i

        if (x+1) % (c_x//50) == 0:
            print('#', end='')
    endTime = datetime.now()
    print('| 100%, done for', endTime - startTime)
    return b


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def binning_c(double[:] pos, long c_x, int P):          # (Chromosome obj, bool for deltas)

    cdef double d_unit = 1.0/P     # unit of 'd' in (0,1)-length chromosome

    cdef long[:] c = np.zeros(P, dtype=int)


    cdef long d_current = 0 # d_current - position in term of d_unit: {0, 1, ... , d_current, ... , P - 1}
    cdef long x
    for x in range(c_x):
        while pos[x] >= (d_current + 1) * d_unit:
            d_current += 1

        c[d_current] += 1                           # for each variant in d_current bin c += 1
    return c



@cython.wraparound(False)
@cython.boundscheck(False)
def FFT(x):
    # Part 2 - FFT

    startTime = datetime.now()
    x = np.array(x)
    if len(x.shape) == 1:
        x = x.reshape(1, -1)

    cdef long p
    p = x.shape[1]
    print('FFT: |', end='')

    cdef long l_x = len(x)
    # # # # # # # # # # # For Output
    step = l_x//50
    power = 1
    odd = 0
    if step == 0:
        step = 1
        power = 50//l_x
        odd = 50%l_x
    # # # # # # # # # # #
    cdef long i
    cdef long j
    cdef long k
    cdef long a
    cdef complex[::1] B = np.zeros(2*p - 1, dtype=complex)
    cdef complex[:,::1] fft_Beauty = np.zeros((2*p - 1, 2*p-1), dtype=complex)
    for i in range(l_x):
        B = fft(np.append(np.array(x[i]), np.zeros(p - 1)))
        for j in range(2*p - 1):
            for k in range(2*p - 1):
                a = j - k
                if a < 0:
                    a = a + 2*p - 1
                fft_Beauty[j][k] = fft_Beauty[j][k] + B[k] * B[j].conjugate() * B[a]  # при замене k на j - будет (ds, d)
    # # # # # # # # # # # # # # # # # For Output
        if i%(step) == 0:
            print('#'*power, end='')
        print('#'*odd, end='')
    # # # # # # # # # # # # # # # # #

    foBeauty = ifft2(fft_Beauty)


    print('| 100%', end='')
    endTime = datetime.now()
    print(', FFT done for', endTime - startTime)

    return foBeauty.real #i don't know why ifft (both 1d and 2d) make it result*n1*n2 (feature of fftw)


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def thld_6(double[:] T, double[:] M, int N, int P):
    cdef double[:,::1] thLD = np.empty((P,P), dtype=np.double)
    cdef int i
    cdef int j
    cdef double d
    cdef double ds
    cdef double l
    for i in range(P):
        for j in range(P):
            l = 1 #1.35 ms  #0.4 rl
            d = i/(l*P) #* length# /(l*P)
            ds = j/(l*P) #* length#/(l*P)
            thLD[i,j] = -(1-M[0]) * (1-M[1]) * exp(-T[1]*(d+ds)) \
    * (  M[1]*(1-M[0])**2 \
         - (2*M[1]**2)*(1-M[0])**2 \
         + M[0]*(1-2*M[0])*exp(-T[0]*(d+ds)) \
         - M[0]*M[1]*(1-M[0]) \
         * ( exp(-T[0]*d) \
             + exp(-T[0]*ds) \
             + (1 - exp(-d) - exp(-ds) \
             + 2*exp(-d-ds))**T[0] )  )
    return thLD
