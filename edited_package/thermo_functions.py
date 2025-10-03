import numpy as np
from .constants import k_B

# Heat capacity from Victor Porée
def partition_func(Eis,T):
    # partition function
    return np.sum(np.exp(-Eis/(k_B*T)))

def Cp_from_CEF(Eis,T):
    def Cp1T(t):
        R = 8.31432  # in J/K per mol
        beta = k_B * t
        Z = partition_func(Eis, t)
        fs = np.sum( (Eis/beta)**2 * np.exp(-Eis/beta) )
        ss = np.sum( (Eis/beta)*np.exp(-Eis/beta) )
        return ((R/Z) * (fs - ss**2/Z))
    return np.array(list(map(Cp1T, T)))
