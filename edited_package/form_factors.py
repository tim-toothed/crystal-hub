"""
Magnetic form factors for neutron scattering calculations.

Implements dipole approximation form factors f(Q) for rare earth ions using
tabulated coefficients from neutron diffraction literature. Form factors
account for spatial distribution of magnetic electrons.

Reference:
    Brown, P.J. in International Tables for Crystallography Vol. C (2006)
"""

import numpy as np
from .constants import ION_NUMS_RARE_EARTH

# Outer Function
def RE_FormFactor(magQ,ion):
    """
    Calculate magnetic form factor for rare earth ion.
    
    Uses dipole approximation: f(Q) = ⟨j₀⟩ + ⟨j₂⟩·C₂
    where C₂ depends on L, S, J quantum numbers.
    
    Args:
        magQ: Momentum transfer magnitude |Q| in inverse Angstroms (scalar)
        ion: Ion symbol (e.g., 'Yb3+', 'Er3+')
    
    Returns:
        float: |f(Q)|² intensity factor
    
    Example:
        >>> Q = 2.5  # Å⁻¹
        >>> ff = RE_FormFactor(Q, 'Yb3+')
        >>> intensity = ff * structure_factor
    
    Note:
        Requires 'RE_formfactors.pck' file with tabulated coefficients
    """
    s = magQ/(4.*np.pi)
    
    coefs0, coefs2 = importRE_FF(ion)
    
    j0 = coefs0[0]*np.exp(-coefs0[1]*s**2) + coefs0[2]*np.exp(-coefs0[3]*s**2)  + \
        coefs0[4]*np.exp(-coefs0[5]*s**2) +coefs0[6]

    j2 = s**2*(coefs2[0]*np.exp(-coefs2[1]*s**2) + coefs2[2]*np.exp(-coefs2[3]*s**2)+\
               coefs2[4]*np.exp(-coefs2[5]*s**2) + coefs2[6])

    S, L, J = ION_NUMS_RARE_EARTH[ion]
     
    j2factor = (J*(J+1.) - S*(S+1.) + L*(L+1.))/(3.*J*(J+1.) + S*(S+1.) - L*(L+1.))
    return (j0 + j2*j2factor)**2

# Inner Function
def importRE_FF(ion):
    """
    Load form factor coefficients from database file.
    
    Reads ⟨j₀⟩ and ⟨j₂⟩ expansion coefficients from RE_formfactors.pck.
    Each is 7-parameter fit: A₁·exp(-a₁·s²) + A₂·exp(-a₂·s²) + A₃·exp(-a₃·s²) + C
    
    Args:
        ion: Ion symbol
    
    Returns:
        tuple: (coefs_j0, coefs_j2) - each is 7-element list [A₁,a₁,A₂,a₂,A₃,a₃,C]
    """
    coefs = [[],[]]
    j=0
    for line in open('RE_formfactors.pck'):
        if not line.startswith(('#', ' ','\n')):
            if line.split(' \t')[0] in ion:
                coefs[j] = [float(i) for i in line.split(' \t')[1:]]
                j+=1
    return coefs[0], coefs[1]