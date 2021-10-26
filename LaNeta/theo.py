import numpy as np


def formula_4_emp(T, M, N, P):
    thLD = np.empty((P,P))

    for i in range(P):
        for j in range(P):

            d = i/P
            ds = j/P

            thLD[i][j] = -M[0] * (1 - M[0]) * (1 - 2*M[0]) \
                       * (1 - 1/(2*N))**T[0] \
                       * (1 - 2/(2*N))**T[0] \
                       * np.exp(-T[0]*(d+ds))
    return thLD



def formula_6(d, ds, T, M, N, P):
    return -(1-M[0]) * (1-M[1]) * np.exp(-T[1]*(d+ds)) \
    * (  M[1]*(1-M[0])**2 \
         - (2*M[1]**2)*(1-M[0])**2 \
         + M[0]*(1-2*M[0])*np.exp(-T[0]*(d+ds)) \
         - M[0]*M[1]*(1-M[0]) \
         * ( np.exp(-T[0]*d) \
             + np.exp(-T[0]*ds) \
             + (1 - np.exp(-d) - np.exp(-ds) \
             + 2*np.exp(-d-ds))**T[0] )  )

rho = 1.6e-9#recombination rate per bp per generation
length=int(1/rho) #1 Morgan chromosome length

def thld_6(T, M, N, P):
    thLD = np.empty((P,P))
    for i in range(P):
        for j in range(P):
            l = 1 #1.35 ms  #0.4 rl
            d = i/(l*P) #* length# /(l*P)
            ds = j/(l*P) #* length#/(l*P)
            thLD[i,j] = formula_6(d, ds, T, M, N, P)
    return thLD



def formula_4(d,ds, D, T, N, drift=False):
    N = N
    T = T[0]+T[1]

    L = np.array(
  [ [ 1,                0,                 0,               0,              0                    ],
    [ 1/(2*N),          (2*N-1)/(2*N),     0,               0,              0                    ],
    [ 1/(2*N),          0,                 (2*N-1)/(2*N),   0,              0                    ],
    [ 1/(2*N),          0,                 0,               (2*N-1)/(2*N),  0                    ],
    [ 1/(4*(N**2)),     (2*N-1)/(2*N),     (2*N-1)/(2*N),   (2*N-1)/(2*N),  (2*N-1)*(2*N-1)/(2*N)] ])
    U = np.array(
  [ [ np.exp(-d-ds),     (1-np.exp(-d))*np.exp(-ds),   (1-np.exp(-ds))*(1-np.exp(-d)),           (1-np.exp(-ds))*np.exp(-d), 0 ],
    [ 0,                 np.exp(-ds),                   0,                                        0,                        (1-np.exp(-ds))],
    [ 0,                 0,                             1-np.exp(-d)-np.exp(-ds)+2*np.exp(-ds-d), 0 ,                        np.exp(-d)+np.exp(-ds)-2*np.exp(-ds-d)],
    [ 0,                 0,                             0,                                        np.exp(-d),                1-np.exp(-d)    ],
    [ 0,                 0,                             0,                                        0,                         1]])


    v = np.array([D[-1][0][0],D[-1][1][1]**2,D[-1][2][2]**2,D[-1][3][3]**2,D[-1][4][4]**3])

    for i in range(T-1,-1,-1):
        v = np.matmul(U, v)
        if drift:
            v = np.matmul(L, v)
        v = np.matmul(D[i], v)

    return np.matmul(np.array([1,-1,-1,-1,2]), v)

def thld_4(D, T, N, P, drift=False):
    scale = P//5-P//200
    print(scale)
    idk = np.zeros((scale,scale))
    for i in range(scale):
        for j in range(scale):
            idk[i][j] = formula_4((i+P//200)/P, (j+P//200)/P, D, T, N, drift)
        if i%30 == 0:
            print('#', end='')
    return idk
