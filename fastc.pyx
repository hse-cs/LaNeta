# cython: language_level=3

cimport cython
from libc.math cimport exp
from libc.math cimport ceil
from libc.math cimport abs
import numpy as np
from datetime import datetime
from numpy.fft import fft, ifft2
from cyvcf2 import VCF
from tqdm import tqdm


from libc.stdio cimport *
from libc.string cimport strchr
from libc.stdlib cimport strtol
from libc.stdlib cimport strtod

# cdef extern from "stdio.h":
#     #FILE * fopen ( const char * filename, const char * mode )
#     FILE *fopen(const char *, const char *)
#     #int fclose ( FILE * stream )
#     int fclose(FILE *)
#     #ssize_t getline(char **lineptr, size_t *n, FILE *stream);
#     ssize_t getline(char **, size_t *, FILE *)

def read_positions(filename, long c_x, long chr_i): #c_x - max count to read, chr_i - chromosome to read

    cdef double[::1] pos = np.empty((c_x), dtype=np.double)

    filename_byte_string = filename.encode("UTF-8")
    cdef char* fname = filename_byte_string

    cdef FILE* cfile
    cfile = fopen(fname, "rb")
    if cfile == NULL:
        raise FileNotFoundError(2, "No such file or directory: '%s'" % filename)

    cdef char * line = NULL
    cdef char * pch = NULL
    cdef size_t l = 0
    cdef ssize_t read
    cdef long ph = 0
    cdef int ws = ord(' ')
    cdef int c = 0

    while True:
        read = getline(&line, &l, cfile)
        if read == -1 or c == c_x: break
        ph = strtol(line, NULL, 10)
        if ph > chr_i: break

        if chr_i == ph:
            pch = strchr(line, ws)
            pch = strchr(pch+1, ws)
            pch = strchr(pch+1, ws)
            if pos[c-1] > strtod(pch, NULL) / 100 : break
            pos[c] = strtod(pch, NULL) / 100
            c += 1


    fclose(cfile)

    return pos[:c], c


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def read_genotype(vcf_gen, unsigned long c_x, int n_min=0, int n_max=0):

    cdef int i = 0
    cdef short[:,::1] y
    cdef unsigned long c = 0
    cdef char[:, ::1] H = np.empty((c_x, n_max-n_min), dtype=np.byte) #but it's not H generally

    for var in tqdm(vcf_gen):
        y = var.genotype.array()[:,:2].copy(order='C')
        for i in range(n_min, n_max):
            if y[i][0] > 1:
                y[i][0] = 1
            if y[i][1] > 1:
                y[i][1] = 1
            H[c][i-n_min] = y[i][0]+y[i][1]

        c = c + 1
        if c == c_x: break
    return H


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def estimate_coef_slow(char[:,::1] H, int n0,
                       char[:,::1] F, int n1,
                       char[:,::1] G, int n2,
                       double[:] pos, int c_x,
                                      double d_unit,
                                      double eps):

    cdef int P_ = int(ceil(pos[c_x-1]/d_unit))

    cdef double[:, ::1] numerator = np.zeros((P_, P_))
    cdef double[:, ::1] denumerator = np.zeros((P_, P_))


    cdef:
        int x = 0, y = 0, z = 0
        int d = 0, ds = 0
        char f = 0, fds = 0
        double f_x = 0, g_x = 0, h_x = 0, sigma_x = 0
        double f_y = 0, g_y = 0, h_y = 0, sigma_y = 0
        double f_z = 0, g_z = 0, h_z = 0, sigma_z = 0
        double Exyzsum = 0

    print(P_)

    for x in range(0, c_x):
        print(pos[x])
        for y in range(x+1, c_x):
            for z in range(y+1, c_x):
                fd = 0
                fds = 0
                for d in range(0, P_):
                    if abs((pos[y] - pos[x]) - d_unit*(d+1)) < eps:
                        fd = 1
                        break
                for ds in range(0, P_):
                    if abs((pos[z]-pos[y]) - d_unit*(ds+1)) < eps:
                        fds = 1
                        break
                if fd == 1 and fds == 1:

                    #Allele freq for H
                    h_x = 0
                    for s_i in range(0, n1):
                        h_x += H[x][s_i]
                    h_x = h_x / n1

                    h_y = 0
                    for s_i in range(0, n1):
                        h_y += H[y][s_i]
                    h_y = h_y / n1

                    h_z = 0
                    for s_i in range(0, n1):
                        h_z += H[z][s_i]
                    h_z = h_z / n1

                    #Expectations sum for H
                    Exyzsum = 0
                    for s_i in range(0, n1):
                        Exyzsum += (H[x][s_i] - h_x) * (H[y][s_i] - h_y) * (H[z][s_i] - h_z)

                    #Deltas for x
                    f_x = 0                                # f_x - impirical allile freq.
                    for s_i in range(0, n1):
                        f_x += F[x][s_i]
                    f_x = f_x / n1

                    g_x = 0                                # f_x - impirical allile freq.
                    for s_i in range(0, n2):
                        g_x += G[x][s_i]
                    g_x = g_x / n2

                    sigma_x = (f_x - g_x)

                    #Deltas for y
                    f_y = 0                                # f_x - impirical allile freq.
                    for s_i in range(0, n1):
                        f_y += F[y][s_i]
                    f_y = f_y / n1

                    g_y = 0                                # f_x - impirical allile freq.
                    for s_i in range(0, n2):
                        g_y += G[y][s_i]
                    g_y = g_y / n2

                    sigma_y = (f_y - g_y)

                    #Deltas for z
                    f_z = 0                                # f_x - impirical allile freq.
                    for s_i in range(0, n1):
                        f_z += F[z][s_i]
                    f_z = f_z / n1

                    g_z = 0                                # f_x - impirical allile freq.
                    for s_i in range(0, n2):
                        g_z += G[z][s_i]
                    g_z = g_z / n2

                    sigma_z = (f_z - g_z)

                    numerator[d][ds] += sigma_x*sigma_y*sigma_z*Exyzsum
                    denumerator[d][ds] += 1
    return numerator, denumerator


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def estimate_affine_term(char[:,:,::1] chroms, int chr_n, int c_x, int n):

    cdef:
        int i=0, j=0, k=0
        int c=0, c_i=0
        double Ex_i=0, Ex_j=0, Ex_k=0
        double at_ijk=0
        double at=0

    for i in range(0, chr_n):
        j = 0
        for j in range(i+1, chr_n):
            if i!=j:
                k = 0
                for k in range(j+1, chr_n):
                    if i!=k:
                        x = 0
                        at_ijk = 0
                        for x in range(0, c_x):
                            s_i = 0
                            Ex_i = 0
                            Ex_j = 0
                            Ex_k = 0
                            for s_i in range(0, n):
                                Ex_i += chroms[i, x, s_i]
                                Ex_j += chroms[j, x, s_i]
                                Ex_k += chroms[k, x, s_i]
                            Ex_i /= n
                            Ex_j /= n
                            Ex_k /= n
                            s_i = 0
                            for s_i in range(0, n):
                                at_ijk += (chroms[i, x, s_i] - Ex_i) * (chroms[j, x, s_i] - Ex_j) * (chroms[k, x, s_i] - Ex_k)
                        at += at_ijk / c_x
    cdef double n0 = n
    at = n0/((n0-1)*(n0-2)) * at
    return at / (chr_n * (chr_n - 1) * (chr_n - 2) / 6)


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def est_deltas(char[:,::1] F, int n1,
                   char[:,::1] G, int n2,
                   double[:] pos, int c_x,):          # (Chromosome obj, bool for deltas)

    cdef int x
    cdef int s_i
    cdef double delta_x = 0
    for x in range(0, c_x):

        F_x = 0                                # f_x - impirical allile freq.
        for s_i in range(0, n1):
            F_x += F[x][s_i]
        F_x = F_x / n1

        G_x = 0                                # f_x - impirical allile freq.
        for s_i in range(0, n2):
            G_x += G[x][s_i]
        G_x = G_x / n2

        delta_x += (F_x - G_x)**2
    delta_x /= c_x
    return delta_x


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def binning_deltas(char[:,::1] F, int n1,
                   char[:,::1] G, int n2,
                   double[:] pos, int c_x,
                   double d_unit):          # (Chromosome obj, bool for deltas)
    # Part 1 - Binning data
    #startTime = datetime.now()
    cdef:
        int bins_n = int(ceil(pos[c_x - 1]/d_unit))
        double[:] delta_b = np.zeros(bins_n, dtype=np.double)
        int d_current = 0 # d_current - position in term of d_unit: {0, 1, ... , d_current, ... , P - 1}
        double F_x = 0
        double G_x = 0
        double sigma_x = 0
        int x = 0, s_i = 0


    for x in range(0, c_x):
        while pos[x] >= (d_current + 1) * d_unit and d_current < bins_n:
            d_current += 1

        if pos[x] <= (d_unit)*bins_n: #may be excessive
            F_x = 0                                # f_x - impirical allile freq.
            for s_i in range(0, n1):
                F_x += F[x][s_i]
            F_x = F_x / n1

            G_x = 0                                # f_x - impirical allile freq.
            for s_i in range(0, n2):
                G_x += G[x][s_i]
            G_x = G_x / n2

            sigma_x = (F_x - G_x) # sigma_x - dif. in allile freq. in two src. pop.

            delta_b[d_current] += sigma_x**2
    #print('Bin deltas: done!')
    return delta_b


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def binning_b_int(char[:,::1] H, int n0,
                  char[:,::1] F, int n1,
                  char[:,::1] G, int n2,
                  double[:] pos, int c_x,
                                 double d_unit):          # (Chromosome obj, bool for deltas)
    # Part 1 - Binning data
    #startTime = datetime.now()
    cdef:
        int d_current = 0 # d_current - position in term of d_unit: {0, 1, ... , d_current, ... , P - 1}
        int bins_n = int(ceil(pos[c_x - 1]/d_unit))

        long f_x, F_x, G_x, sigma_x
        int x, s_i

        long[:,::1] b = np.zeros((n0, bins_n), dtype=np.int_)

    for x in range(0, c_x):
        while pos[x] >= (d_current + 1) * d_unit and d_current < bins_n:
            d_current += 1

        f_x = 0                                # f_x - impirical allile freq.
        for s_i in range(0, n0):
            f_x += H[x][s_i]


        F_x = 0                                # f_x - impirical allile freq.
        for s_i in range(0, n1):
            F_x += F[x][s_i]


        G_x = 0                                # f_x - impirical allile freq.
        for s_i in range(0, n2):
            G_x += G[x][s_i]

        sigma_x = (n2*F_x - n1*G_x) # sigma_x - dif. in allile freq. in two src. pop.

        for s_i in range(0, n0):
            b[s_i, d_current] += sigma_x * (H[x, s_i]*n0 - f_x)      # computing b_i[d_current] for every i

    #endTime = datetime.now()
    #print('b bin done for', endTime - startTime)
    return b


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def binning_b(char[:,::1] H, int n0,
              char[:,::1] F, int n1,
              char[:,::1] G, int n2,
              double[:] pos, int c_x,
                             double d_unit):
    #startTime = datetime.now()
    cdef:
        int bins_n = int(ceil(pos[c_x - 1]/d_unit))
        cdef int d_current = 0 # d_current - position in term of d_unit: {0, 1, ... , d_current, ... , P - 1}

        double f_x, F_x, G_x, sigma_x
        int x, s_i

        cdef double[:,::1] b = np.zeros((n0, bins_n), dtype=np.double)

    for x in range(0, c_x):
        while pos[x] >= (d_current + 1) * d_unit and d_current < bins_n:
            d_current += 1

        if pos[x] <= (d_unit)*bins_n: #may be excessive

            f_x = 0                                # f_x - impirical allile freq.
            for s_i in range(0, n0):
                f_x += H[x][s_i]
            f_x = f_x / n0

            F_x = 0                                # f_x - impirical allile freq.
            for s_i in range(0, n1):
                F_x += F[x][s_i]
            F_x = F_x / n1

            G_x = 0                                # f_x - impirical allile freq.
            for s_i in range(0, n2):
                G_x += G[x][s_i]
            G_x = G_x / n2

            sigma_x = (F_x - G_x) # sigma_x - dif. in allile freq. in two src. pop.

            for s_i in range(0, n0):
                b[s_i, d_current] += sigma_x * (H[x, s_i] - f_x)      # computing b_i[d_current] for every i

    #endTime = datetime.now()
    #print('b bin done for', endTime - startTime)
    return b


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def binning_c(double[:] pos, long c_x, double d_unit):          # (Chromosome obj, bool for deltas)

    cdef int bins_n = int(ceil(pos[c_x - 1]/d_unit))
    cdef long[:] c = np.zeros(bins_n, dtype=int)
    cdef long d_current = 0 # d_current - position in term of d_unit: {0, 1, ... , d_current, ... , P - 1}
    cdef long x
    for x in range(0, c_x):
        while pos[x] >= (d_current + 1) * d_unit and d_current < bins_n:
            d_current += 1

        c[d_current] += 1                           # for each variant in d_current bin c += 1
    return c


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def noFFT(double[:,::1] A, int s, int P): #improve without P

    cdef:
        double[:,::1] Anp = np.zeros((s, 2*P - 1)), B = np.zeros((P, P))
        int d = 0, ds = 0, w = 0,i = 0

    Anp[:, :P] = A

    for i in range(0, s):
        for d in range(0, P):
            for ds in range(0, P - d):
                for w in range(0, P):
                    B[d, ds] += Anp[i, w]*Anp[i, w+d]*Anp[i, w+d+ds]
    print('noftend')
    return B



@cython.wraparound(False)
@cython.boundscheck(False)
def FFT(x):

    x = np.array(x)
    if len(x.shape) == 1:
        x = x.reshape(1, -1)

    cdef:
        long p = x.shape[1], i, j, k, a, l_x = len(x)
        complex[::1] B = np.zeros(2*p - 1, dtype=complex)
        complex[:,::1] fft_Beauty = np.zeros((2*p - 1, 2*p-1), dtype=complex)

    for i in range(0, l_x):
        B = fft(np.append(np.array(x[i]), np.zeros(p - 1)))
        for j in range(2*p - 1):
            for k in range(2*p - 1):
                a = j - k
                if a < 0:
                    a = a + 2*p - 1
                fft_Beauty[j][k] = fft_Beauty[j][k] + B[k] * B[j].conjugate() * B[a]  # при замене k на j - будет (ds, d)
    foBeauty = ifft2(fft_Beauty)


    return foBeauty.real


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def thld_6(double[:] T, double[:] M, int N, double d_unit, int P_):

    cdef:
        double[:,::1] thLD = np.empty((P_,P_), dtype=np.double)
        int i, j
        double d, ds, l

    for i in range(0, P_):
        for j in range(0, P_):
            #max_l = 1 #1.35 ms  #0.4 rl
            d = (d_unit*i) #* length# /(l*P)
            ds = (d_unit*j) #* length#/(l*P)
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


def formula_4(d,ds, D, T, N, drift=False):
    T = T[0]+T[1]

    L = np.array(
  [ [ (4*(N**2)),     0,                 0,               0,              0                    ],
    [ (2*N),          (2*N-1)*(2*N),     0,               0,              0                    ],
    [ (2*N),          0,                 (2*N-1)*(2*N),   0,              0                    ],
    [ (2*N),          0,                 0,               (2*N-1)*(2*N),  0                    ],
    [ 1,              (2*N-1),           (2*N-1),         (2*N-1),        (2*N-1)*(2*N-1)      ] ])
    U = np.array(
  [ [ np.exp(-d-ds),     (1-np.exp(-d))*np.exp(-ds),   (1-np.exp(-ds))*(1-np.exp(-d)),           (1-np.exp(-ds))*np.exp(-d), 0 ],
    [ 0,                 np.exp(-ds),                   0,                                        0,                        (1-np.exp(-ds))],
    [ 0,                 0,                             1-np.exp(-d)-np.exp(-ds)+2*np.exp(-ds-d), 0 ,                        np.exp(-d)+np.exp(-ds)-2*np.exp(-ds-d)],
    [ 0,                 0,                             0,                                        np.exp(-d),                1-np.exp(-d)    ],
    [ 0,                 0,                             0,                                        0,                         1]])


    v = np.array([D[-1][0][0],D[-1][1][1],D[-1][2][2],D[-1][3][3],D[-1][4][4]])

    for i in range(T-1,-1,-1):
        if drift:
            A = np.matmul(D[i], L)
            A = np.matmul(A, U)
        else:
            A = np.matmul(D[i], U)
        v = np.matmul(A, v)
        if drift:
            v /= (4*N**2)

    return np.matmul(np.array([1,-1,-1,-1,2]), v)


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def thld_4(D, T, int N, double d_unit, int P_, drift=False):
    cdef:
        double[:,::1] thLD = np.empty((P_,P_), dtype=np.double)
        int i, j
    for i in range(0, P_):
        for j in range(0, P_):
            thLD[i][j] = formula_4(d_unit * i, d_unit * j, D, T, N, drift)
    return thLD
