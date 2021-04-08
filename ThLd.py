import algorithm
import gen
import estimator
import numpy as np

class ThLd:
    def __init__(self, data=''):
        if data == 'ms':
            params = np.loadtxt('tap.txt')
            self.TrueTs = params[0]
            self.TrueMs = params[1]
            self.data = gen.twopulse(self.TrueTs, self.TrueMs)
        else:
            self.data = []

    def alg(self, deltas=True, coef=True):
        if len(self.data)!=0:
            c, d = algorithm.alg(self.data, deltas=deltas, coef=coef)
        else:
            print('No data')
        if deltas == True:
            self.deltas = d
        if coef == True:
            self.coef = c

    def load(self, deltas='', coef='', data=''):
        if deltas != '':
            self.deltas = np.load(deltas)
            print(deltas, '- succesfully loaded!')
        if coef != '':
            self.coef = np.load(coef)
            print(coef, '- succesfully loaded!')
        #if data != '':
        #    self.data = np.load(data)
        #    print(data, '- succesfully loaded!')

    def estim(self):
        print('estimating parameteres ...')
        res = estimator.estimate(self.deltas, self.coef)
        print(res)
        self.Ts = res.x
        print(self.Ts)
