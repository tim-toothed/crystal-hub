from typing import Union, Dict

import numpy as np


#######################
# Free-Ion State Data
#######################

### These values are taken from the free ion states at 
### https://physics.nist.gov/PhysRefData/Elements/per_noframes.html

# S,L,J = quantum numbers of Ion
# s = spin quantum number
# l = orbital angular momentum quantum number
# j = total angular momentum quantum number

# Transition metals with [S, L] values
ION_NUMS_TRANS_METAL = {
    'Ag2+': [1/2, 2], 'Ag3+': [1, 3],
    'Cd3+': [1/2, 2],
    'Co2+': [3/2, 3], 'Co3+': [2, 2], 'Co6+': [3/2, 3],
    'Cr2+': [2, 2], 'Cr3+': [3/2, 3], 'Cr4+': [1, 3], 'Cr5+': [1/2, 2],
    'Cu2+': [1/2, 2.],
    'Fe2+': [2, 2], 'Fe3+': [5/2, 0],
    'Hf2+': [1, 3], 'Hf3+': [1/2, 2],
    'Mn2+': [5/2, 0], 'Mn3+': [2, 2], 'Mn4+': [3/2, 3], 'Mn5+': [1, 3], 'Mn6+': [1/2, 2],
    'Mo2+': [2, 2], 'Mo3+': [3/2, 3], 'Mo4+': [1, 3], 'Mo5+': [1/2, 2],
    'Nb3+': [1, 3],
    'Ni2+': [1., 3.], 'Ni3+': [3/2, 3.],
    'Pd2+': [1, 3], 'Pd3+': [3/2, 3], 'Pd4+': [2, 2],
    'Re3+': [2, 2], 'Re4+': [3/2, 3], 'Re6+': [1/2, 2],
    'Rh2+': [3/2, 3], 'Rh3+': [2, 2], 'Rh4+': [5/2, 0],
    'Ru2+': [2, 2], 'Ru3+': [5/2, 0], 'Ru4+': [2, 2], 'Ru6+': [1, 3],
    'Ta2+': [3/2, 3], 'Ta3+': [1, 3], 'Ta4+': [1/2, 2],
    'Tc4+': [3/2, 3],
    'Ti2+': [1, 3], 'Ti3+': [1/2, 2],
    'V2+': [3/2, 3], 'V3+': [1, 3], 'V4+': [1/2, 2],
    'W2+': [2, 2], 'W3+': [3/2, 3], 'W4+': [1, 3], 'W5+': [1/2, 2], 'W6+': [0, 1],
    'Y2+': [1/2, 2],
    'Zr+': [3/2, 3], 'Zr2+': [1, 3], 'Zr3+': [1/2, 2],
}

# Rare Earths with [S, L, J] values
ION_NUMS_RARE_EARTH = {
    'Ce3+': [0.5, 3., 2.5], 'Pr3+': [1., 5., 4.], 'Nd3+': [1.5, 6., 4.5],
    'Pm3+': [2., 6., 4.], 'Sm3+': [2.5, 5, 2.5], 'Eu3+': [3, 3, 0],
    'Gd3+': [7/2, 0, 7/2],
    'Tb3+': [3., 3., 6.], 'Dy3+': [2.5, 5., 7.5], 'Ho3+': [2., 6., 8.],
    'Er3+': [1.5, 6., 7.5], 'Tm3+': [1., 5., 6.], 'Yb3+': [0.5, 3., 3.5],
    'U4+': [1., 5., 4.], 'U3+': [1.5, 6., 4.5],
}


#######################
# Constants for converting between Wybourne and Stevens Operators
#######################

WYBOURNE_STEVENS_CONSTS = {
    2: {0: 1/2, 1: np.sqrt(6), 2: np.sqrt(6)/2},
    4: {0: 1/8, 1: np.sqrt(5)/2, 2: np.sqrt(10)/4, 3: np.sqrt(35)/2, 4: np.sqrt(70)/8},
    6: {0: 1/16, 1: np.sqrt(42)/8, 2: np.sqrt(105)/16, 3: np.sqrt(105)/8,
         4: 3*np.sqrt(14)/16, 5: 3*np.sqrt(77)/8, 6: np.sqrt(231)/16}
}

######################
# Tesseral Constants
######################

# Define constants as a nested dictionary for clean access
TESSERAL_CONSTANTS: Dict[int, Dict[int, float]] = {
    0: {0: 0.282094791773878},
    
    1: {-1: 0.48860251190292, 0: 0.48860251190292, 1: 0.48860251190292},
    
    2: {-2: 0.54627421529604, -1: 1.09254843059208, 0: 0.31539156525252, 
        1: 1.09254843059208, 2: 0.54627421529604},
    
    3: {-3: 0.590043589926644, -2: 1.44530572132028, -1: 0.457045799464466,
        0: 0.373176332590115, 1: 0.457045799464466, 2: 1.44530572132028, 
        3: 0.590043589926644},
    
    4: {-4: 0.625835735449176, -3: 1.77013076977993, -2: 0.47308734787878,
        -1: 0.669046543557289, 0: 0.105785546915204, 1: 0.669046543557289,
        2: 0.47308734787878, 3: 1.77013076977993, 4: 0.625835735449176},
    
    5: {-5: 0.65638205684017, -4: 2.07566231488104, -3: 0.489238299435251,
        -2: 2.39676839248666, -1: 0.452946651195697, 0: 0.116950322453424,
        1: 0.452946651195697, 2: 2.39676839248666, 3: 0.489238299435251,
        4: 2.07566231488104, 5: 0.65638205684017},
    
    6: {-6: 0.683184105191914, -5: 2.36661916223175, -4: 0.504564900728724,
        -3: 0.921205259514924, -2: 0.460602629757462, -1: 0.582621362518732,
        0: 0.0635692022676284, 1: 0.582621362518732, 2: 0.460602629757462,
        3: 0.921205259514924, 4: 0.504564900728724, 5: 2.36661916223175,
        6: 0.683184105191914}
}

# <- ХУЙНЯ ФУНКЦИЯ - УБРАТЬ, ПЕРЕДЕЛАТЬ, СЖЕЧЬ
def Constant(n: int, m: int) -> float:
    """Returns the constant in front of the tesseral harmonic."""
    return TESSERAL_CONSTANTS[n][m]

# Function for "Calculate tesseral harmonic function values."
def calculate_tesseral_harmonic(n: int, m: int, x: Union[float, np.ndarray], 
             y: Union[float, np.ndarray], z: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Calculate tesseral harmonic function values.
    These functions have been cross-checked with Mathematica.
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    r_sq = r**2
    r_cubed = r**3
    r_4 = r**4
    r_6 = r**6

    def _tesseral_dispatch(n: int, m: int, x: float, y: float, z: float, 
                           r: float, r_sq: float, r_cubed: float, r_4: float, r_6: float) -> float:
        """Dispatch to the appropriate tesseral harmonic calculation."""
        if n == 0 and m == 0:
            return 1.0
        
        elif n == 1:
            if m == 1: return x / r
            elif m == 0: return z / r
            elif m == -1: return y / r
        
        elif n == 2:
            if m == -2: return 2 * x * y / r_sq
            elif m == -1: return y * z / r_sq
            elif m == 0: return (3 * z**2 - r_sq) / r_sq
            elif m == 1: return x * z / r_sq
            elif m == 2: return (x**2 - y**2) / r_sq
        
        elif n == 3:
            if m == -3: return (3 * x**2 * y - y**3) / r_cubed
            elif m == -2: return (2 * x * y * z) / r_cubed
            elif m == -1: return y * (5 * z**2 - r_sq) / r_cubed
            elif m == 0: return z * (5 * z**2 - 3 * r_sq) / r_cubed
            elif m == 1: return x * (5 * z**2 - r_sq) / r_cubed
            elif m == 2: return z * (x**2 - y**2) / r_cubed
            elif m == 3: return (x**3 - 3 * x * y**2) / r_cubed
        
        elif n == 4:
            if m == -4: return 4 * (x**3 * y - x * y**3) / r_4
            elif m == -3: return (3 * x**2 * y - y**3) * z / r_4
            elif m == -2: return 2 * x * y * (7 * z**2 - r_sq) / r_4
            elif m == -1: return y * z * (7 * z**2 - 3 * r_sq) / r_4
            elif m == 0: return (35 * z**4 - 30 * z**2 * r_sq + 3 * r_4) / r_4
            elif m == 1: return x * z * (7 * z**2 - 3 * r_sq) / r_4
            elif m == 2: return (x**2 - y**2) * (7 * z**2 - r_sq) / r_4
            elif m == 3: return (x**3 - 3 * x * y**2) * z / r_4
            elif m == 4: return (x**4 - 6 * x**2 * y**2 + y**4) / r_4
        
        elif n == 5:
            return 0.0  # Skipping n=5 as per original comment
        
        elif n == 6:
            if m == -6: return (6 * x**5 * y - 20 * x**3 * y**3 + 6 * x * y**5) / r_6
            elif m == -5: return (5 * x**4 * y - 10 * x**2 * y**3 + y**5) * z / r_6
            elif m == -4: return 4 * (x**3 * y - x * y**3) * (11 * z**2 - r_sq) / r_6
            elif m == -3: return (3 * x**2 * y - y**3) * (11 * z**3 - 3 * z * r_sq) / r_6
            elif m == -2: return 2 * x * y * (33 * z**4 - 18 * z**2 * r_sq + r_4) / r_6
            elif m == -1: return y * z * (33 * z**4 - 30 * z**2 * r_sq + 5 * r_4) / r_6
            elif m == 0: return (231 * z**6 - 315 * z**4 * r_sq + 105 * z**2 * r_4 - 5 * r**6) / r_6
            elif m == 1: return x * z * (33 * z**4 - 30 * z**2 * r_sq + 5 * r_4) / r_6
            elif m == 2: return (x**2 - y**2) * (33 * z**4 - 18 * z**2 * r_sq + r_4) / r_6
            elif m == 3: return (x**3 - 3 * x * y**2) * (11 * z**3 - 3 * z * r_sq) / r_6
            elif m == 4: return (x**4 - 6 * x**2 * y**2 + y**4) * (11 * z**2 - r_sq) / r_6
            elif m == 5: return (x**5 - 10 * x**3 * y**2 + 5 * x * y**4) * z / r_6
            elif m == 6: return (x**6 - 15. * x**4 * y**2 + 15. * x**2 * y**4 - y**6) / r_6
        
        raise ValueError(f"Tesseral harmonic not implemented for n={n}, m={m}")
    
    # Use a dispatch pattern based on (n, m) tuple
    value = _tesseral_dispatch(n, m, x, y, z, r, r_sq, r_cubed, r_4, r_6) 
    
    from . import TESSERAL_CONSTANTS # <- короче так будет, я хз как по-умному

    return TESSERAL_CONSTANTS[n][m] * value


######################
# Multiplicative factor from Hutchings, Table VII
######################

# This functions are not used in any package files <- по идее это можно нахуй убрать, либо наоборот убрать дикты ниже

def PFalpha(L: float, S: float, l: float, halffilled: bool = True) -> float:
    """
    Calculate the PFalpha multiplicative factor from Hutchings, Table VII.
    
    Args:
        L: Orbital angular momentum quantum number
        S: Spin angular momentum quantum number  
        l: Orbital quantum number
        halffilled: Whether the shell is half-filled (default: True)
    
    Returns:
        PFalpha factor
    """
    numerator = 2.0 * (2.0 * l + 1.0 - 4.0 * S)
    denominator = (2.0 * l - 1.0) * (2.0 * l + 3.0) * (2.0 * L - 1.0)
    
    factor = numerator / denominator
    return factor if halffilled else -factor


def PFbeta(L: float, S: float, l: float, halffilled: bool = True) -> float:
    """
    Calculate the PFbeta multiplicative factor from Hutchings, Table VII.
    
    Args:
        L: Orbital angular momentum quantum number
        S: Spin angular momentum quantum number
        l: Orbital quantum number  
        halffilled: Whether the shell is half-filled (default: True)
    
    Returns:
        PFbeta factor
    """
    # Break down the complex numerator into manageable parts
    spin_term = 2.0 * l + 1.0 - 4.0 * S
    bracket_term = (-7.0 * (l - 2.0 * S) * (l - 2.0 * S + 1.0) + 
                    3.0 * (l - 1.0) * (l + 2.0))
    
    numerator = 3.0 * spin_term * bracket_term
    
    # Break down the complex denominator
    l_terms = ((2.0 * l - 3.0) * (2.0 * l - 1.0) * 
               (2.0 * l + 3.0) * (2.0 * l + 5.0))
    L_terms = (L - 1.0) * (2.0 * L - 1.0) * (2.0 * L - 3.0)
    
    denominator = l_terms * L_terms
    
    factor = numerator / denominator
    return factor if halffilled else -factor


def PFgamma(L: float, nvalence: int) -> float:
    """
    Calculate the PFgamma multiplicative factor for rare earth ions.
    
    Note: We assume l=6 because only l=6 shells have a gamma term.
    
    Args:
        L: Orbital angular momentum quantum number for the multiplet
        nvalence: Number of valence electrons in the f-shell
    
    Returns:
        PFgamma factor
    """
    # Calculate O_06 expectation value for the multiplet
    X_multiplet = L * (L + 1.0)
    multiplet_expectation = (
        231.0 * L**6 
        - (315.0 * X_multiplet - 735.0) * L**4 
        + (105.0 * X_multiplet**2 - 525.0 * X_multiplet + 294.0) * L**2 
        - 5.0 * X_multiplet**3 + 40.0 * X_multiplet**2 - 60.0 * X_multiplet
    )
    
    # Gamma6 constant from integration over spherical harmonics in l=3, m=3 state
    GAMMA6_CONSTANT = -4.0 / 3861.0
    
    # Calculate individual electron wave function expectation values
    # For f-electrons (l=3), possible m_l values range from -3 to +3
    magnetic_quantum_numbers = np.tile(np.arange(-3, 4), 2)
    
    individual_electron_sum = 0.0
    for i in range(nvalence):
        Lz = magnetic_quantum_numbers[i]
        X_electron = 3.0 * (3.0 + 1.0)  # l=3 for f-electrons
        
        electron_contribution = (
            231.0 * Lz**6 
            - (315.0 * X_electron - 735.0) * Lz**4 
            + (105.0 * X_electron**2 - 525.0 * X_electron + 294.0) * Lz**2 
            - 5.0 * X_electron**3 + 40.0 * X_electron**2 - 60.0 * X_electron
        )
        individual_electron_sum += electron_contribution
    
    factor = (individual_electron_sum / multiplet_expectation) * GAMMA6_CONSTANT
    
    return factor

########################################################################
# The following lookup table was generated from the functions above.
# This was done to save time in computation steps.

def LStheta(ion,n): # <- не дает ли неточностей эти консты с 13 знаками после запятой??
    LSThet = {}
    LSThet['Sm3+'] = [0.0148148148148, 0.0003848003848, -2.46666913334e-05]
    LSThet['Pm3+'] = [0.0040404040404, 0.000122436486073, 1.12121324243e-05]
    LSThet['Nd3+'] = [-0.0040404040404, -0.000122436486073, -1.12121324243e-05]
    LSThet['Ce3+'] = [-0.0444444444444, 0.0040404040404, -0.001036001036]
    LSThet['Dy3+'] = [-0.0148148148148, -0.0003848003848, 2.46666913334e-05]
    LSThet['Ho3+'] = [-0.0040404040404, -0.000122436486073, -1.12121324243e-05]
    LSThet['Tm3+'] = [0.0148148148148, 0.0003848003848, -2.46666913334e-05]
    LSThet['Pr3+'] = [-0.0148148148148, -0.0003848003848, 2.46666913334e-05]
    LSThet['U4+'] = [-0.0148148148148, -0.0003848003848, 2.46666913334e-05] #same as Pr3+
    LSThet['Er3+'] = [0.0040404040404, 0.000122436486073, 1.12121324243e-05]
    LSThet['Tb3+'] = [-0.0444444444444, 0.0040404040404, -0.001036001036]
    LSThet['Yb3+'] = [0.0444444444444, -0.0040404040404, 0.001036001036]
    if isinstance(ion, str):
        return LSThet[ion][int(n/2-1)] # <- почему нельзя было сразу посчитать?


# Multiplicative factor for rare earth ground state multiplet
# from Hutchings, Table VI
 # Cross-checked by hand with calculator.

def theta(ion,n): # <- здесь я вообще хз искренне
    Thet = {}
    Thet['Ce3+'] = [-2./(5*7), 2./(3*3*5*7), 0]
    Thet['Pr3+'] = [-2.*2*13/(3*3*5*5*11), -2.*2/(3*3*5*11*11), 2.**4*17/(3**4*5*7*11**2*13)]
    Thet['U4+'] = [-2.*2*13/(3*3*5*5*11), -2.*2/(3*3*5*11*11), 2.**4*17/(3**4*5*7*11**2*13)] #same as Pr3+
    Thet['Nd3+'] = [-7./(3**2*11**2) , -2.**3*17/(3**3*11**3*13), -5.*17*19/(3**3*7*11**3*13**2)]
    Thet['U3+'] = [-7./(3**2*11**2) , -2.**3*17/(3**3*11**3*13), -5.*17*19/(3**3*7*11**3*13**2)] #same as Nd3+
    Thet['Pm3+'] = [2*7./(3*5*11**2), 2.**3*7*17/(3**3*5*11**3*13), 2.**3*17*19/(3**3*7*11**2*13**2)]
    Thet['Sm3+'] = [13./(3**2*5*7) , 2.*13/(3**3*5*7*11), 0]
    Thet['Tb3+'] = [-1./(3**2*11), 2./(3**3*5*11**2), -1./(3**4*7*11**2*13)]
    Thet['Dy3+'] = [-2./(3**2*5*7) , -2.**3/(3**3*5*7*11*13), 2.*2/(3**3*7*11**2*13**2)]
    Thet['Ho3+'] = [-1./(2*3*3*5*5), -1./(2*3*5*7*11*13), -5./(3**3*7*11**2*13**2)]
    Thet['Er3+'] = [2.*2/(3*3*5*5*7) , 2./(3.**2*5*7*11*13), 2.*2*2/(3.**3*7*11**2*13**2)]
    Thet['Tm3+'] = [1./(3**2*11) , 2.**3/(3**4*5*11**2), -5./(3**4*7*11**2*13)]
    Thet['Yb3+'] = [2./(3**2*7) , -2./(3*5*7*11), 2.*2/(3**3*7*11*13)]
    return Thet[ion][int(n/2-1)]


#######################
# Spin Orbit Coupling Constants
#######################

# Constants in cm^-1 (converted to meV when used)
# Sources:
# - https://pubs.acs.org/doi/abs/10.1021/acs.jpca.8b09218 (EXPT or PMT preferred)
# - Fe3+, Mn2+, and Ru3+ from https://arxiv.org/pdf/cond-mat/0505214.pdf
# - W6+ from https://doi.org/10.1016/0166-1280(95)04297-0

# <- тоже как бы нахуй не надо, вообще можно оставить только SPIN_ORBIT_COUPLING_CONSTANTS, опять вычисления показывает
SPIN_ORBIT_COUPLING_CM = {
    'Ti2+': 121,
    'Ti3+': 154,
    'V2+': 168,
    'V3+': 201,
    'V4+': 248,
    'Cr2+': 234,
    'Cr3+': 276,
    'Cr4+': 329,
    'Cr5+': 383,
    'Mn2+': 322,
    'Mn3+': 358,
    'Mn4+': 405,
    'Mn5+': 479,
    'Mn6+': 542,
    'Fe2+': -401,
    'Fe3+': 476,
    'Co2+': -527,
    'Co3+': -619,
    'Co6+': 787,
    'Ni2+': -643,
    'Ni3+': -749,
    'Cu2+': -829,
    'Y2+': 290,
    'Zr+': 347,
    'Zr2+': 428,
    'Zr3+': 500,
    'Nb3+': 677,
    'Mo2+': 714,
    'Mo3+': 836,
    'Mo4+': 972,
    'Mo5+': 1031,
    'Tc4+': 1150,
    'Ru2+': -942,
    'Ru3+': 1177,
    'Ru4+': 1350,
    'Ru6+': 1700,
    'Rh2+': -1194,
    'Rh3+': -1360,
    'Rh4+': 1350,
    'Pd2+': -1293,
    'Pd3+': -1211,
    'Pd4+': -1830,
    'Ag2+': -1843,
    'Ag3+': -1930,
    'Cd3+': -2325,
    'Hf2+': 1370,
    'Hf3+': 1877,
    'Ta2+': 2000,
    'Ta3+': 2254,
    'Ta4+': 2643,
    'W2+': 2500,
    'W3+': 2686,
    'W4+': 3108,
    'W5+': 3483,
    'W6+': 5070,
    'Re3+': 5800,
    'Re4+': 3300,
    'Re6+': 4398,
    'Os2+': -2200,
    'Os4+': 4000,
    'Os5+': 4500,
    'Os6+': 5200,
    'Os7+': 5200,
    'Ir2+': -3900,
    'Ir3+': -3050,
    'Ir5+': 5500,
    'Ir6+': 6000,
    'Pt2+': -4500,
    'Pt4+': -3700,
    'Au3+': -5200,
}

#HC_CONVERSION = 1.23984193e-1  # meV*cm conversion factor

# Convert to meV and create the final dictionary
#SPIN_ORBIT_COUPLING_CONSTANTS = {ion: coupling_cm * HC_CONVERSION 
#             for ion, coupling_cm in SPIN_ORBIT_COUPLING_CM.items()}

SPIN_ORBIT_COUPLING_CONSTANTS = {
    'Ti2+': 15.002087353,
    'Ti3+': 19.093565722,
    'V2+': 20.829344424000002,
    'V3+': 24.920822793000003,
    'V4+': 30.748079864,
    'Cr2+': 29.012301162,
    'Cr3+': 34.219637268,
    'Cr4+': 40.790799497,
    'Cr5+': 47.485945919,
    'Mn2+': 39.922910146,
    'Mn3+': 44.386341094,
    'Mn4+': 50.213598165,
    'Mn5+': 59.388428447,
    'Mn6+': 67.199432606,
    'Fe2+': -49.717661393,
    'Fe3+': 59.016475868,
    'Co2+': -65.339669711,
    'Co3+': -76.746215467,
    'Co6+': 97.57555989100001,
    'Ni2+': -79.721836099,
    'Ni3+': -92.864160557,
    'Cu2+': -102.78289599700001,
    'Y2+': 35.955415970000004,
    'Zr+': 43.022514971,
    'Zr2+': 53.065234604000004,
    'Zr3+': 61.9920965,
    'Nb3+': 83.937298661,
    'Mo2+': 88.52471380200001,
    'Mo3+': 103.650785348,
    'Mo4+': 120.51263559600001,
    'Mo5+': 127.82770298300001,
    'Tc4+': 142.58182195,
    'Ru2+': -116.793109806,
    'Ru3+': 145.929395161,
    'Ru4+': 167.37866055,
    'Ru6+': 210.7731281,
    'Rh2+': -148.03712644200002,
    'Rh3+': -168.61850248000002,
    'Rh4+': 167.37866055,
    'Pd2+': -160.311561549,
    'Pd3+': -150.144857723,
    'Pd4+': -226.89107319000001,
    'Ag2+': -228.502867699,
    'Ag3+': -239.28949249000001,
    'Cd3+': -288.26324872500004,
    'Hf2+': 169.85834441,
    'Hf3+': 232.718330261,
    'Ta2+': 247.968386,
    'Ta3+': 279.460371022,
    'Ta4+': 327.69022209900004,
    'W2+': 309.9604825,
    'W3+': 333.021542398,
    'W4+': 385.342871844,
    'W5+': 431.83694421900003,
    'W6+': 628.59985851,
    'Re3+': 719.1083194,
    'Re4+': 409.1478369,
    'Re6+': 545.282480814,
    'Os2+': -272.7652246,
    'Os4+': 495.936772,
    'Os5+': 557.9288685,
    'Os6+': 644.7178036,
    'Os7+': 644.7178036,
    'Ir2+': -483.5383527,
    'Ir3+': -378.15178865,
    'Ir5+': 681.9130615,
    'Ir6+': 743.905158,
    'Pt2+': -557.9288685,
    'Pt4+': -458.7415141,
    'Au3+': -644.7178036
    }

#######################
# Radial Integrals
#######################

#Import Radial Integrals
RADIAL_INTEGRALS_RARE_EARTH = {
    'Ce3+': [0.51, 0.0132, -0.0294, 1.456, 5.437, 42.26],
    'Dy3+': [0.527, -0.0199, -0.0316, 0.849, 1.977, 10.44],
    'Er3+': [0.544, -0.0427, -0.031, 0.773, 1.677, 8.431],
    'Eu3+': [0.52, 0.0033, -0.0319, 0.997, 2.638, 15.34],
    'Gd3+': [0.521, -0.0031, -0.0318, 0.942, 2.381, 13.36],
    'Ho3+': [0.534, -0.0306, -0.0313, 0.81, 1.816, 9.345],
    'Lu3+': [0.588, -0.0902, -0.0294, 0.682, 1.353, 6.441],
    'Nd3+': [0.518, 0.013, -0.031, 1.222, 3.875, 26.12],
    'Pm3+': [0.519, 0.0109, -0.0314, 1.135, 3.366, 21.46],
    'Pr3+': [0.515, 0.0138, -0.0301, 1.327, 4.537, 32.65],
    'Sm3+': [0.519, 0.0077, -0.0317, 1.061, 2.964, 17.99],
    'Tb3+': [0.523, -0.0107, -0.0318, 0.893, 2.163, 11.75],
    'Tm3+': [0.554, -0.0567, -0.0306, 0.74, 1.555, 7.659],
    'Yb3+': [0.571, -0.0725, -0.03, 0.71, 1.448, 7.003]
    }

RADIAL_INTEGRALS_TRANS_METAL = {
    'Ag2+': [0.55, 0.559], 'Ag3+': [0.499, 0.437],
    'Au3+': [0.635, 0.663],
    'Cd3+': [0.732, 1.739],
    'Co2+': [0.36, 0.31], 'Co3+': [0.302, 0.2], 'Co6+': [0.206, 0.081],
    'Cr2+': [0.515, 0.605], 'Cr3+': [0.416, 0.362], 'Cr4+': [0.347, 0.238], 'Cr5+': [0.302, 0.172],
    'Cu2+': [0.295, 0.21],
    'Fe2+': [0.402, 0.38], 'Fe3+': [0.334, 0.24],
    'Hf2+': [1.372, 3.284], 'Hf3+': [1.162, 2.212],
    'Ir2+': [0.811, 1.135], 'Ir3+': [0.736, 0.89], 'Ir5+': [0.632, 0.617], 'Ir6+': [0.594, 0.533],
    'Mn2+': [0.452, 0.475], 'Mn3+': [0.371, 0.293], 'Mn4+': [0.316, 0.201], 'Mn5+': [0.274, 0.144], 'Mn6+': [0.244, 0.11],
    'Mo2+': [0.949, 1.633], 'Mo3+': [0.825, 1.162], 'Mo4+': [0.735, 0.887], 'Mo5+': [0.667, 0.708],
    'Nb3+': [0.935, 1.482],
    'Ni2+': [0.325, 0.256], 'Ni3+': [0.275, 0.168],
    'Os2+': [0.885, 1.357], 'Os4+': [0.731, 0.845], 'Os5+': [0.677, 0.707], 'Os6+': [0.634, 0.606],
    'Pd2+': [0.617, 0.705], 'Pd3+': [0.555, 0.54], 'Pd4+': [0.545, 0.519],
    'Pt2+': [0.747, 0.961], 'Pt4+': [0.633, 0.635],
    'Re3+': [0.868, 1.237], 'Re4+': [0.789, 0.986], 'Re6+': [0.678, 0.692],
    'Rh2+': [0.668, 0.819], 'Rh3+': [0.61, 0.649], 'Rh4+': [0.546, 0.497],
    'Ru2+': [0.744, 1.012], 'Ru3+': [0.674, 0.79], 'Ru4+': [0.599, 0.596], 'Ru6+': [0.51, 0.409],
    'Ta2+': [1.209, 2.545], 'Ta3+': [1.047, 1.8], 'Ta4+': [0.934, 1.375],
    'Tc4+': [0.674, 0.752],
    'Ti2+': [0.699, 1.075], 'Ti3+': [0.538, 0.587],
    'V2+': [0.595, 0.793], 'V3+': [0.47, 0.456], 'V4+': [0.388, 0.292],
    'W2+': [1.08, 2.025], 'W3+': [0.95, 1.482], 'W4+': [0.857, 1.16], 'W5+': [0.785, 0.946], 'W6+': [0.113, 0.027],
    'Y2+': [1.535, 4.173],
    'Zr+': [4.102, 25.077], 'Zr2+': [1.282, 2.94], 'Zr3+': [1.072, 1.932]
}

def calculate_radial_integral_RE(ion,n):
    """Returns the radial integral of a rare earth ion plus self-shielding.
    Comes out in units of Bohr radius"""

    from . import RADIAL_INTEGRALS_RARE_EARTH # <- короче так будет, я хз как по-умному

    if ion == 'U4+':
        U4r = {2: 2.042, 4: 7.632, 6: 47.774}  # from Freeman, Desclaux, Lander, and Faber, PRB (1976), Table I
        return U4r[n]
    elif ion == 'U3+':
        U3r = {2: 2.346, 4: 10.906, 6: 90.544}  # from Freeman, Desclaux, Lander, and Faber, PRB (1976), Table I
        return U3r[n]
    else:
        shielding = 1- RADIAL_INTEGRALS_RARE_EARTH[ion][int(n/2-1)]
        return RADIAL_INTEGRALS_RARE_EARTH[ion][int(n/2-1) + 3] * shielding

def calculate_radial_integral_TM(ion,n):

    from . import RADIAL_INTEGRALS_TRANS_METAL  # <- короче так будет, я хз как по-умному

    """Returns the radial integral of a transition ion.
    The listed constants are in AA, so we convert to Bohr Radii"""
    BohrRadius= 0.5291772109 # Official NIST value in units of /AA
    return RADIAL_INTEGRALS_TRANS_METAL[ion][int(n/2-1)]/(BohrRadius**n)


#######################
# Half-filled Ions
#######################

ION_HALF_FILLED = ['Mn2+',
            'Fe2+','Fe3+',
            'Co2+','Co3+',
            'Ni2+','Ni3+',
            'Cu2+',
            'Ru2+',
            'Rh2+','Rh3+',
            'Ag2+','Ag3+',
            'Cd3+',
            'Os2+',
            'Ir2+','Ir3+',
            'Au3+']

ION_NOT_HALF_FILLED = ['Ti2+','Ti3+',
                'V2+','V3+','V4+',
                'Cr2+','Cr3+','Cr4+','Cr5+',
                'Mn3+','Mn4+','Mn5+','Mn6+',
                'Co6+',
                'Y2+',
                'Zr2+','Zr3+',
                'Nb3+',
                'Mo2+','Mo3+','Mo4+','Mo5+',
                'Tc4+',
                'Ru3+','Ru4+','Ru6+',
                'Rh4+',
                'Pd2+','Pd3+','Pd4+',
                'Hf2+','Hf3+',
                'Ta2+','Ta3+','Ta4+',
                'W2+','W3+','W4+','W5+','W6+',
                'Re3+','Re4+','Re6+',
                'Os4+','Os5+','Os6+','Os7+',
                'Ir5+','Ir6+',
                'Pt2+','Pt4+']

def is_half_filled(ion):
    '''determine if given ion has a half-filled shell or not.'''

    from . import ION_HALF_FILLED, ION_NOT_HALF_FILLED  # <- короче так будет, я хз как по-умному

    if ion in ION_HALF_FILLED:
        return True
    elif ion in ION_NOT_HALF_FILLED:
        return False
    else:
        raise ValueError('{} is not a known ion for PyCrystalField.'.format(ion))

#######################
# Other Constants
#######################

ahc = 1.43996e4  #Constant to get the energy in units of meV = alpha*hbar*c
a0 = 0.52917721067    #Bohr radius in \AA
muB = 5.7883818012e-2  # meV/T
k_B = 8.6173303e-2  # meV/K