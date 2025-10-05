import numpy as np
from .constants import k_B

# used outside the package

def partition_func(Eis,T):
    """
    Calculates the partition function Z(T) for a quantum system

    Args:
        Eis: Energy eigenvalues (array) in meV
        T: Temperature (float or array) in Kelvin

    Returns:
        Partition function value: Z = Σ exp(-Eᵢ/kᵦT)
    """

    return np.sum(np.exp(-Eis/(k_B*T)))

def Cp_from_CEF(Eis,T):
    """
    Calculates molar heat capacity from crystal field energy levels

    Args:
        Eis: Energy eigenvalues (array) in meV
        T: Temperature array in Kelvin

    Returns:
        Heat capacity array Cₚ(T) in J/(K·mol)
    """
    
    def Cp1T(t):
        R = 8.31432  # in J/K per mol
        beta = k_B * t
        Z = partition_func(Eis, t)
        fs = np.sum( (Eis/beta)**2 * np.exp(-Eis/beta) )
        ss = np.sum( (Eis/beta)*np.exp(-Eis/beta) )
        return ((R/Z) * (fs - ss**2/Z))
    return np.array(list(map(Cp1T, T)))
