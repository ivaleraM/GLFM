import numpy as np

class Hidden:

    def __init__(self, Z, B, Theta, mu, wu, s2y):
        self.Z = Z
        self.B = B
        self.Theta = Theta
        self.mu = mu
        self.wu = wu
        self.s2y = s2y
