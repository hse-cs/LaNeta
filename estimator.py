from scipy.optimize import minimize
from algorithm import P
import numpy as np
import theo
import gen

M1 = 0.2
M2 = 0.2
T1 = 1
T2 = 1

def metr(var, D, a_coef):
    print(var)
    params = gen.Params(var[0], var[1], M1, M2)
    th_E = theo.thld_6(params)
    return np.sqrt((((D[:P//5,:P//5])*th_E[:P//5,:P//5] - a_coef[:P//5,:P//5])**2).sum())

def estimate(deltas, a_coef):
    var = np.array([T1, T2])
    res = minimize(metr, var, (deltas, a_coef), method='BFGS')
    return res
