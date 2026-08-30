"""
Parameters for He*-He ULRRM potential calculation.
"""

# IMPORTS
import numpy as np

class Params:

    def __init__(self):

        self.n = 10
        self.n_min = self.n - 2
        self.n_max = self.n + 2
        self.n_arr = np.arange(self.n_min, self.n_max)

        self.l = 8
        self.l_min = 3

        self.basis = [
            (n, l)
            for n in self.n_arr
            for l in range(self.l_min, n)
        ]

        self.m = 0

        self.R_arr = np.linspace(180,190,100) # Ground He perturber position relative to He* ionic core
        #self.R_arr = [190]

        self.phi = 0
        self.theta = 0

        self.classically_forbidden = False