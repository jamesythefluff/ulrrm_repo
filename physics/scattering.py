import numpy as np

class Scattering:

    def __init__(self, scattering_data):

        self.data = scattering_data

    def a_s(self, k):

        # Ignoring classically forbidden region - improve later
        if k == None:
            return 0
        
        a_s = - np.tan(self.data.calculate_delta_0(k)) / k

        return a_s

    def a_p(self, k):

        # Ignoring classically forbidden region - improve later
        if k == None:
            return 0

        a_p = - np.tan(self.data.calculate_delta_1(k)) / k**3

        return a_p