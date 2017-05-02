import numpy as np

class Params:

    def __init__(self, missing=-1, alpha=1, bias=0, s2u=0.01, \
            s2B=1, Niter=1000, maxK=50, verbose=1, t=[], t_1=[], dt_1=[]):
        self.missing = missing
        self.alpha = alpha
        self.bias = bias
        self.s2u = s2u
        self.s2B = s2B
        self.Niter = Niter
        self.maxK = maxK
        self.verbose = verbose
        self.t = t
        self.t_1 = t_1
        self.dt_1 = dt_1
