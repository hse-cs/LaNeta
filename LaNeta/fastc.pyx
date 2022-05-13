# cython: language_level=3

cimport cython
from libc.math cimport exp
from libc.math cimport ceil
from libc.math cimport abs
import numpy as np
from datetime import datetime
from numpy.fft import fft, ifft2
from scipy.stats import binom
from cyvcf2 import VCF


from libc.stdio cimport *
from libc.string cimport strchr
from libc.string cimport strcmp
from libc.string cimport strtok
from libc.stdlib cimport strtol
from libc.stdlib cimport strtod

# cdef extern from "stdio.h":
#     #FILE * fopen ( const char * filename, const char * mode )
#     FILE *fopen(const char *, const char *)
#     #int fclose ( FILE * stream )
#     int fclose(FILE *)
#     #ssize_t getline(char **lineptr, size_t *n, FILE *stream);
#     ssize_t getline(char **, size_t *, FILE *)


def read_positions(filename, long c_x, str chr_name): #c_x - max count to read, chr_i - chromosome to read

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
    cdef char *  ph = NULL
    cdef int ws = ord(' ')
    cdef int c = 0

    while True:
        read = getline(&line, &l, cfile)
        if read == -1 or c == c_x: break
        #ph = strtol(line, NULL, 10)
        ph = strtok(line, ' ')
        #if ph > chr_i: break
        if ph == chr_name.encode('UTF-8'):
            pch = strtok(NULL, ' ')
            pch = strtok(NULL, ' ')
            pch = strtok(NULL, ' ')
            #if pos[c-1] > strtod(pch, NULL) / 100 : break
            #print(pch)
            pos[c] = strtod(pch, NULL) / 100
            c += 1


    fclose(cfile)

    return pos[:c], c


@cython.wraparound(False)
@cython.boundscheck(False)
def read_vcf_pos(vcf_gen, long c_x): #c_x - max count to read, chr_i - chromosome to read

    cdef double[::1] pos = np.empty((c_x), dtype=np.double)
    cdef long j = 0
    for var in vcf_gen:
        pos[j] = float(var.INFO.get('pos'))
        j += 1

    return pos


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def read_vcf_genotype(vcf_gen, long c_x, long[:] H_samples, long[:] F_samples, long[:] G_samples, long[:] inds = None):

    cdef int i = 0
    cdef short[:,::1] y
    cdef long c = 0
    cdef int n0 = len(H_samples) , n1 = len(F_samples), n2 = len(G_samples)
    cdef char[:, ::1] H = np.empty((c_x, n0), dtype=np.byte) #but it's not H generally
    cdef char[:, ::1] F = np.empty((c_x, n1), dtype=np.byte) #but it's not F generally
    cdef char[:, ::1] G = np.empty((c_x, n2), dtype=np.byte) #but it's not G generally

    cdef long id = 0

    cdef long j = 0

    for var in vcf_gen:
        if j == c_x: break
        if inds == None or c == inds[j]:
            y = var.genotype.array()[:,:2].copy(order='C')
            for i in range(n0):
                id = H_samples[i]
                if y[id][0] > 1:
                    y[id][0] = 1
                if y[id][1] > 1:
                    y[id][1] = 1
                H[j][i] = y[id][0]+y[id][1]
            for i in range(n1):
                id = F_samples[i]
                if y[id][0] > 1:
                    y[id][0] = 1
                if y[id][1] > 1:
                    y[id][1] = 1
                F[j][i] = y[id][0]+y[id][1]
            for i in range(n2):
                id = G_samples[i]
                if y[id][0] > 1:
                    y[id][0] = 1
                if y[id][1] > 1:
                    y[id][1] = 1
                G[j][i] = y[id][0]+y[id][1]
            j += 1
        c = c + 1
        #if c == c_x: break
    return H, F, G


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
cdef double estimate_freq(char[:] H, int n):
    cdef double freq = 0
    cdef int i = 0
    cdef int freq_n = 0
    for i in range(n):
        if H[i] >= 0:
            freq+=H[i]
            freq_n+=1
    if freq_n!=0:
        freq /= freq_n
    return freq


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
cdef double estimate_freq_double(double[:] H, int n):
    cdef double freq = 0
    cdef int i = 0
    cdef int freq_n = 0
    for i in range(n):
        if H[i] >= 0:
            freq+=H[i]
            freq_n+=1
    if freq_n!=0:
        freq /= freq_n
    return freq


cdef missing_freq(double h, double f, double m):
    cdef double freq = (h - f*m)/(1-m)
    if freq<0:
        freq = 0
    if freq>1:
        freq = 1
    return freq

def missing_del(double h, double f, double m):
    cdef double freq = (f - h)/(m)
    return freq



@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def estimate_coef_slow(double[:,::1] H, int n0,
                       double[:,::1] F, int n1,
                       double[:,::1] G, int n2,
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
def estimate_affine_term_v2(double[:,:,::1] H_chroms, double[:,:,::1] F_chroms, double[:,:,::1] G_chroms, int chr_n, int c_x,
                            int n0, int n1, int n2,
                            char missing_ref, double m):

    cdef:
        int i=0, j=0, k=0
        int c=0, c_i=0
        double hx_i=0, hx_j=0, hx_k=0
        double fx_i=0, fx_j=0, fx_k=0
        double gx_i=0, gx_j=0, gx_k=0
        double Deltax_i=0, Deltax_j=0, Deltax_k=0
        double at_ijk=0, at_ijk_s=0
        double at=0
        int n_ = 0
        int x = 0

    rng = np.random.default_rng()

    for x in range(0, c_x):
        i, j, k = rng.choice(chr_n, size=3, replace=False)

        hx_i = estimate_freq_double(H_chroms[i, x], n0)
        hx_j = estimate_freq_double(H_chroms[j, x], n0)
        hx_k = estimate_freq_double(H_chroms[k, x], n0)
        fx_i = estimate_freq_double(F_chroms[i, x], n1)
        gx_i = estimate_freq_double(G_chroms[i, x], n2)
        fx_j = estimate_freq_double(F_chroms[j, x], n1)
        gx_j = estimate_freq_double(G_chroms[j, x], n2)
        fx_k = estimate_freq_double(F_chroms[k, x], n1)
        gx_k = estimate_freq_double(G_chroms[k, x], n2)
        #print(h_x, f_x, g_x)
        if missing_ref==1:
            fx_i = (hx_i*n0 + fx_i*n1) / (n0+n1)
            fx_i=missing_freq(fx_i, gx_i, m)
            fx_j = (hx_j*n0 + fx_j*n1) / (n0+n1)
            fx_j=missing_freq(fx_j, gx_j, m)
            fx_k = (hx_k*n0 + fx_k*n1) / (n0+n1)
            fx_k=missing_freq(fx_k, gx_k, m)
        if missing_ref==2:
            gx_i = (hx_i*n0 + gx_i*n2) / (n0+n2)
            gx_i=missing_freq(gx_i, fx_i, 1-m)
            gx_j = (hx_j*n0 + gx_j*n2) / (n0+n2)
            gx_j=missing_freq(gx_j, fx_j, 1-m)
            gx_k = (hx_k*n0 + gx_k*n2) / (n0+n2)
            gx_k=missing_freq(gx_k, fx_k, 1-m)
        Deltax_i = fx_i-gx_i
        Deltax_j = fx_j-gx_j
        Deltax_k = fx_k-gx_k
        s_i = 0
        n_ = 0
        at_ijk_s = 0
        for s_i in range(0, n0):
            if H_chroms[i, x, s_i] >= 0 and H_chroms[i, x, s_i] >= 0 and H_chroms[i, x, s_i] >= 0 and fx_i+hx_i>0.1 and fx_j+hx_j>0.1 and fx_k+hx_k>0.1:
                at_ijk_s += (H_chroms[i, x, s_i] - hx_i) * (H_chroms[j, x, s_i] - hx_j) * (H_chroms[k, x, s_i] - hx_k) * \
                            Deltax_i * Deltax_j * Deltax_k
                n_ += 1
        if n_ != 0:
            at_ijk_s = (at_ijk_s / float(n_)) * float(n0)
        at = (at*x + at_ijk_s)/(x+1)
        print(at)
    at = n0 / ((n0-1.0)*(n0-2.0)) * at
    print(at)
    return at


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def estimate_affine_term(double[:,:,::1] H_chroms, double[:,:,::1] F_chroms, double[:,:,::1] G_chroms, int chr_n, int c_x,
                         int n0, int n1, int n2,
                         char missing_ref, double m):

    cdef:
        int i=0, j=0, k=0
        int c=0, c_i=0
        double hx_i=0, hx_j=0, hx_k=0
        double fx_i=0, fx_j=0, fx_k=0
        double gx_i=0, gx_j=0, gx_k=0
        double Deltax_i=0, Deltax_j=0, Deltax_k=0
        double at_ijk=0, at_ijk_s=0
        double at=0
        int n_ = 0

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
                            hx_i = estimate_freq_double(H_chroms[i, x], n0)
                            hx_j = estimate_freq_double(H_chroms[j, x], n0)
                            hx_k = estimate_freq_double(H_chroms[k, x], n0)
                            fx_i = estimate_freq_double(F_chroms[i, x], n1)
                            gx_i = estimate_freq_double(G_chroms[i, x], n2)
                            fx_j = estimate_freq_double(F_chroms[j, x], n1)
                            gx_j = estimate_freq_double(G_chroms[j, x], n2)
                            fx_k = estimate_freq_double(F_chroms[k, x], n1)
                            gx_k = estimate_freq_double(G_chroms[k, x], n2)
                            #print(h_x, f_x, g_x)
                            if missing_ref==1:
                                fx_i = (hx_i*n0 + fx_i*n1) / (n0+n1)
                                fx_i=missing_freq(fx_i, gx_i, m)
                                fx_j = (hx_j*n0 + fx_j*n1) / (n0+n1)
                                fx_j=missing_freq(fx_j, gx_j, m)
                                fx_k = (hx_k*n0 + fx_k*n1) / (n0+n1)
                                fx_k=missing_freq(fx_k, gx_k, m)
                            if missing_ref==2:
                                gx_i = (hx_i*n0 + gx_i*n2) / (n0+n2)
                                gx_i=missing_freq(gx_i, fx_i, 1-m)
                                gx_j = (hx_j*n0 + gx_j*n2) / (n0+n2)
                                gx_j=missing_freq(gx_j, fx_j, 1-m)
                                gx_k = (hx_k*n0 + gx_k*n2) / (n0+n2)
                                gx_k=missing_freq(gx_k, fx_k, 1-m)
                            Deltax_i = fx_i-gx_i
                            Deltax_j = fx_j-gx_j
                            Deltax_k = fx_k-gx_k
                            s_i = 0
                            n_ = 0
                            at_ijk_s = 0
                            for s_i in range(0, n0):
                                if H_chroms[i, x, s_i] >= 0 and H_chroms[i, x, s_i] >= 0 and H_chroms[i, x, s_i] >= 0:
                                    at_ijk_s += (H_chroms[i, x, s_i] - hx_i) * (H_chroms[j, x, s_i] - hx_j) * (H_chroms[k, x, s_i] - hx_k) * \
                                                Deltax_i * Deltax_j * Deltax_k
                                    n_ += 1
                            if n_ != 0:
                                at_ijk += (at_ijk_s / float(n_)) * float(n0)
                        at += at_ijk / c_x
                        print(at)
        print(i)
    print(at)
    at = n0 / ((n0-1.0)*(n0-2.0)) * at
    print(at)
    return at / (chr_n * (chr_n - 1.0) * (chr_n - 2.0) / 6)


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def est_deltas(char[:,::1] H, int n0,
               char[:,::1] F, int n1,
               char[:,::1] G, int n2,
               double freq_filter, int c_x,
               double m=0, missing_ref=0):          # (Chromosome obj, bool for deltas)

    cdef int x = 0
    cdef double delta_x = 0
    cdef int c = 0

    cdef int norm_c = 0
    cdef double ph_x = 0
    cdef double f_x = 0, h_x = 0, g_x = 0

    for x in range(0, c_x):
        f_x = estimate_freq(F[x], n1)
        h_x = estimate_freq(H[x], n0)
        g_x = estimate_freq(G[x], n2)
        if missing_ref==1:
            #f_x = (h_x*n0 + f_x*n1) / (n0+n1)
            f_x=missing_freq(f_x, g_x, m)
        if missing_ref==2:
            #g_x = (h_x*n0 + g_x*n2) / (n0+n2)
            g_x=missing_freq(g_x, f_x, 1-m)

        if g_x+f_x >= freq_filter:
            delta_x += (f_x - g_x)**2
            norm_c += 1
    delta_x /= norm_c
    #print(delta_x)
    return delta_x


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def binning_deltas(char[:,::1] H, int n0,
                   char[:,::1] F, int n1,
                   char[:,::1] G, int n2,
                   double[:] pos, int c_x,
                   double d_unit,
                   double m=0,
                   char missing_ref=0):          # (Chromosome obj, bool for deltas)
    # Part 1 - Binning data
    #startTime = datetime.now()
    cdef:
        int bins_n = int(ceil(pos[c_x - 1]/d_unit))
        double[:] delta_b = np.zeros(bins_n, dtype=np.double)
        int d_current = 0 # d_current - position in term of d_unit: {0, 1, ... , d_current, ... , P - 1}
        double  h_x = 0, f_x = 0, g_x = 0, delta_x = 0
        int x = 0, s_i = 0


    for x in range(0, c_x):
        while pos[x] >= (d_current + 1) * d_unit and d_current < bins_n:
            d_current += 1

        if pos[x] <= (d_unit)*bins_n: #may be excessive
            h_x = estimate_freq(H[x], n0)
            f_x = estimate_freq(F[x], n1)
            g_x = estimate_freq(G[x], n2)

            if missing_ref==1:
                f_x=missing_freq(f_x, g_x, m)
            if missing_ref==2:
                g_x=missing_freq(g_x, f_x, 1-m)

            delta_x = (f_x - g_x) # sigma_x - dif. in allile freq. in two src. pop.

            delta_b[d_current] += delta_x**2
    #print('Bin deltas: done!')
    return delta_b






# def get_center(A, d_unit):
#     m = ceil(A[-1]/d_unit)
#     C = np.arange(m) * d_unit + (d_unit/2)
#     return C
#
#
#
# def get_bin(dummy, d_unit, r):
#     C = get_center(dummy, d_unit)
#     print(C)
#
#     c1 = 0
#     c2 = 1
#
#     R = np.zeros(len(C))
#
#     for x in dummy:
#         if c1 != -1:
#             if abs(x-C[c1]) <= r:
#                 R[c1] += 1
#             elif x-C[c1] > r:
#                 c1 = c1 + 2
#
#         if c2 != -1:
#             if abs(x-C[c2]) <= r:
#                 R[c2] += 1
#             elif x-C[c2] > r:
#                 c2 = c2 + 2
#
#         if c2 >= len(C):
#             c2 = -1
#         if c1 >= len(C):
#             c1 = -1
#     return R


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def binning_b_double(double[:,::1] H, int n0,
                     double[:,::1] F, int n1,
                     double[:,::1] G, int n2,
                     long[:] c,
                     double[:] pos, int c_x,
                     double d_unit, double r,
                     double m=0,
                     char missing_ref=0):          # (Chromosome obj, bool for deltas)
    # Part 1 - Binning data
    #startTime = datetime.now()
    cdef:
        int prob = 0

        int d_current = 0 # d_current - position in term of d_unit: {0, 1, ... , d_current, ... , P - 1}
        int bins_n = int(ceil(pos[c_x - 1]/d_unit))

        double h_x, f_x, g_x, delta_x
        int x, s_i

        double[:,::1] b = np.zeros((n0, bins_n), dtype=np.double)

        double[:,::1] C_new = np.zeros((n0, bins_n), dtype=np.double)

    cdef double[:] C = np.arange(bins_n, dtype=np.double) * d_unit + (d_unit/2)

    cdef int cc = 0
    cdef int rc = 0
    cdef int lc = 0

    for x in range(0, c_x):

        h_x = estimate_freq_double(H[x], n0)
        f_x = estimate_freq_double(F[x], n1)
        g_x = estimate_freq_double(G[x], n2)
        #print(h_x, f_x, g_x)
        if missing_ref==1:
            f_x=missing_freq(f_x, g_x, m)
        if missing_ref==2:
            g_x=missing_freq(g_x, f_x, 1-m)
            if g_x == 1 or g_x == 0:
                prob+=1
        delta_x = (f_x - g_x)

        while pos[x] >= (d_unit)*(cc+1):
            cc += 1

        rc = 1
        lc = -1

        for s_i in range(0, n0):
            if H[x, s_i] >= 0 and h_x+f_x>0.1: #
                b[s_i, cc] += delta_x * (H[x, s_i] - h_x)      # computing b_i[d_current] for every i
                C_new[s_i, cc] += 1

        while cc+rc < bins_n and -(pos[x] - C[cc+rc]) <= r:
            for s_i in range(0, n0):
                if H[x, s_i] >= 0 and h_x+f_x>0.1:
                    b[s_i, cc+rc] += delta_x * (H[x, s_i] - h_x)      # computing b_i[d_current] for every i
                    C_new[s_i, cc+rc] += 1
            rc += 1

        while cc+lc >= 0 and pos[x] - C[cc+lc] <= r:
            for s_i in range(0, n0):
                if H[x, s_i] >= 0 and h_x+f_x>0.1:
                    b[s_i, cc+lc] += delta_x * (H[x, s_i] - h_x)      # computing b_i[d_current] for every i
                    C_new[s_i, cc+lc] += 1
            lc -= 1


    for s_i in range(0, n0):
        for d_current in range(0, bins_n):
            if C_new[s_i, d_current]!=0:
                b[s_i, d_current] /= C_new[s_i, d_current]
                b[s_i, d_current] *= float(c[d_current])
    #endTime = datetime.now()
    #print('b bin done for', endTime - startTime)
    print('bad missing freq:',prob)
    return b


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def binning_b_int(char[:,::1] H, int n0,
                  char[:,::1] F, int n1,
                  char[:,::1] G, int n2,
                  double freq_filter, long[:] c,
                  double[:] pos, int c_x,
                  double d_unit, double r,
                  double m=0,
                  char missing_ref=0):          # (Chromosome obj, bool for deltas)
    # Part 1 - Binning data
    #startTime = datetime.now()
    cdef:
        int prob = 0

        int d_current = 0 # d_current - position in term of d_unit: {0, 1, ... , d_current, ... , P - 1}
        int bins_n = int(ceil(pos[c_x - 1]/d_unit))

        double h_x, f_x, g_x, delta_x
        int x, s_i

        double[:,::1] b = np.zeros((n0, bins_n), dtype=np.double)

        double[:,::1] C_new = np.zeros((n0, bins_n), dtype=np.double)

    cdef double[:] C = np.arange(bins_n, dtype=np.double) * d_unit + (d_unit/2)

    cdef int cc = 0
    cdef int rc = 0
    cdef int lc = 0

    for x in range(0, c_x):

        h_x = estimate_freq(H[x], n0)
        f_x = estimate_freq(F[x], n1)
        g_x = estimate_freq(G[x], n2)
        #print(h_x, f_x, g_x)
        if missing_ref==1:
            f_x=missing_freq(f_x, g_x, m)
        if missing_ref==2:
            g_x=missing_freq(g_x, f_x, 1-m)
            if g_x == 1 or g_x == 0:
                prob+=1
        delta_x = (f_x - g_x)

        while pos[x] >= (d_unit)*(cc+1):
            cc += 1

        rc = 1
        lc = -1

        for s_i in range(0, n0):
            if H[x, s_i] >= 0 and f_x+g_x >= freq_filter: #
                b[s_i, cc] += delta_x * (H[x, s_i] - h_x)      # computing b_i[d_current] for every i
                C_new[s_i, cc] += 1

        while cc+rc < bins_n and -(pos[x] - C[cc+rc]) <= r:
            for s_i in range(0, n0):
                if H[x, s_i] >= 0 and f_x+g_x >= freq_filter:
                    b[s_i, cc+rc] += delta_x * (H[x, s_i] - h_x)      # computing b_i[d_current] for every i
                    C_new[s_i, cc+rc] += 1
            rc += 1

        while cc+lc >= 0 and pos[x] - C[cc+lc] <= r:
            for s_i in range(0, n0):
                if H[x, s_i] >= 0 and f_x+g_x >= freq_filter:
                    b[s_i, cc+lc] += delta_x * (H[x, s_i] - h_x)      # computing b_i[d_current] for every i
                    C_new[s_i, cc+lc] += 1
            lc -= 1


    for s_i in range(0, n0):
        for d_current in range(0, bins_n):
            if C_new[s_i, d_current]!=0:
                b[s_i, d_current] /= C_new[s_i, d_current]
                b[s_i, d_current] *= float(c[d_current])
    #endTime = datetime.now()
    #print('b bin done for', endTime - startTime)
    print('bad missing freq:',prob)
    return b

#
# @cython.cdivision(True)
# @cython.wraparound(False)
# @cython.boundscheck(False)
# def binning_b_int_int(char[:,::1] H, int n0,
#                   char[:,::1] F, int n1,
#                   char[:,::1] G, int n2,
#                   long[:] c,
#                   double[:] pos, int c_x,
#                   double d_unit, double r,
#                   double m=0,
#                   char missing_ref=0):          # (Chromosome obj, bool for deltas)
#     # Part 1 - Binning data
#     #startTime = datetime.now()
#     cdef:
#         int prob = 0
#
#         int d_current = 0 # d_current - position in term of d_unit: {0, 1, ... , d_current, ... , P - 1}
#         int bins_n = int(ceil(pos[c_x - 1]/d_unit))
#
#         double h_x, f_x, g_x, delta_x
#         int x, s_i
#
#         double[:,::1] b = np.zeros((n0, bins_n), dtype=np.double)
#
#         double[:,::1] C_new = np.zeros((n0, bins_n), dtype=np.double)
#
#     cdef double[:] C = np.arange(bins_n, dtype=np.double) * d_unit + (d_unit/2)
#
#     cdef int cc = 0
#     cdef int rc = 0
#     cdef int lc = 0
#
#     for x in range(0, c_x):
#
#         h_x = estimate_freq(H[x], n0)
#         f_x = estimate_freq(F[x], n1)
#         g_x = estimate_freq(G[x], n2)
#         #print(h_x, f_x, g_x)
#         if missing_ref==1:
#             f_x=missing_freq(h_x, g_x, m)#missing_freq_ex(int(h_x*n0), n0, int(g_x*n2), n2, hst_h, hst_g, m)
#         if missing_ref==2:
#             g_x=missing_freq(h_x, f_x, 1-m)#missing_freq_ex(int(h_x*n0), n0, int(f_x*n1), n1, hst_h, hst_f, 1-m)
#             if g_x > 1 or g_x < 0:
#                 prob+=1
#         delta_x = (f_x - g_x)
#
#         while pos[x] >= (d_unit)*(cc+1):
#             cc += 1
#
#         rc = 1
#         lc = -1
#
#         for s_i in range(0, n0):
#             if H[x, s_i] >= 0: #and f_x>0.05 and h_x>0.05:
#                 b[s_i, cc] += delta_x * (H[x, s_i] - h_x)      # computing b_i[d_current] for every i
#                 C_new[s_i, cc] += 1
#
#         while cc+rc < bins_n and -(pos[x] - C[cc+rc]) <= r:
#             for s_i in range(0, n0):
#                 if H[x, s_i] >= 0: #and f_x>0.05 and h_x>0.05:
#                     b[s_i, cc+rc] += delta_x * (H[x, s_i] - h_x)      # computing b_i[d_current] for every i
#                     C_new[s_i, cc+rc] += 1
#             rc += 1
#
#         while cc+lc >= 0 and pos[x] - C[cc+lc] <= r:
#             for s_i in range(0, n0):
#                 if H[x, s_i] >= 0: #and f_x>0.05 and h_x>0.05:
#                     b[s_i, cc+lc] += delta_x * (H[x, s_i] - h_x)      # computing b_i[d_current] for every i
#                     C_new[s_i, cc+lc] += 1
#             lc -= 1
#
#
#     for s_i in range(0, n0):
#         for d_current in range(0, bins_n):
#             if C_new[s_i, d_current]!=0:
#                 b[s_i, d_current] /= C_new[s_i, d_current]
#                 b[s_i, d_current] *= float(c[d_current])
#     #endTime = datetime.now()
#     #print('b bin done for', endTime - startTime)
#     return b


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


# @cython.cdivision(True)
# @cython.wraparound(False)
# @cython.boundscheck(False)
def binning_c(double[:] pos, long c_x, double d_unit, double r):          # (Chromosome obj, bool for deltas)

    cdef int bins_n = int(ceil(pos[c_x - 1]/d_unit))
    cdef long[:] c = np.zeros(bins_n, dtype=int)
    cdef long d_current = 0 # d_current - position in term of d_unit: {0, 1, ... , d_current, ... , P - 1}
    cdef long x

    cdef double[:] C = np.arange(bins_n, dtype=np.double) * d_unit + (d_unit/2)

    cdef int cc = 0

    for x in range(0, c_x):
        while pos[x] >= (d_unit)*(cc+1):
            cc += 1

        rc = 1
        lc = -1

        c[cc] += 1

        while cc+rc < bins_n and -(pos[x] - C[cc+rc]) <= r:
            c[cc+rc] += 1
            rc += 1

        while cc+lc >= 0 and pos[x] - C[cc+lc] <= r:
            c[cc+lc] += 1
            lc -= 1
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
    X = np.array(x)
    for i in range(0, l_x):
        B = fft(np.append(X[i], np.zeros(p - 1)))
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
            A = np.matmul(L, D[i])
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
