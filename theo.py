import numpy as np
from algorithm import P


def formula_4_emp(params):
    thLD = np.empty((P,P))

    for i in range(P):
        for j in range(P):

            d = i/P
            ds = j/P

            thLD[i][j] = -params.M1 * (1 - params.M1) * (1 - 2*params.M1) \
                       * (1 - 1/(2*params.N))**params.T \
                       * (1 - 2/(2*params.N))**params.T \
                       * np.exp(-params.T*(d+ds))
    return thLD



def formula_6(d, ds, params):
    return -(1-params.M1) * (1-params.M2) * np.exp(-params.T2*(d+ds)) \
    * (  params.M2*(1-params.M1)**2 \
         - (2*params.M2**2)*(1-params.M1)**2 \
         + params.M1*(1-2*params.M1)*np.exp(-params.T1*(d+ds)) \
         - params.M1*params.M2*(1-params.M1) \
         * ( np.exp(-params.T1*d) \
             + np.exp(-params.T1*ds) \
             + (1 - np.exp(-d) - np.exp(-ds) \
             + 2*np.exp(-d-ds))**params.T1 )  )

def thld_6(params):
    thLD = np.empty((P,P))
    for i in range(P):
        for j in range(P):
            d = i/P
            ds = j/P
            thLD[i,j] = formula_6(d, ds, params)
    return thLD



def formula_4(d,ds, params, drift=False):
    m1 = params.M1
    m2 = params.M2
    N = params.N
    T1 = params.T1
    T2 = params.T2
    T = params.T

    D1 = np.diag((1-m1,(1-m1)**2,(1-m1)**2,(1-m1)**2,(1-m1)**3))
    D2 = np.diag((1-m2,(1-m2)**2,(1-m2)**2,(1-m2)**2,(1-m2)**3))
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


    v = np.array([1-m1,(1-m1)**2,(1-m1)**2,(1-m1)**2,(1-m1)**3])

    for i in range(T-1,-1,-1):
        v = np.matmul(U, v)
        if drift:
            v = np.matmul(L, v)
        if i == T2:
            v = np.matmul(D2, v)

    return np.matmul(np.array([1,-1,-1,-1,2]), v)

def thld_4(params, drift=False):
    idk = np.empty((P,P))
    for i in range(P):
        for j in range(P-i):
            idk[i][j] = formula_4(i/P, j/P, params, drift)
        if i%30 == 0:
            print('#', end='')
    return idk
