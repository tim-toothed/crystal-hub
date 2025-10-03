import numpy as np
from .constants import ION_NUMS_RARE_EARTH

def importRE_FF(ion):
    coefs = [[],[]]
    j=0
    for line in open('RE_formfactors.pck'):
        if not line.startswith(('#', ' ','\n')):
            if line.split(' \t')[0] in ion:
                coefs[j] = [float(i) for i in line.split(' \t')[1:]]
                j+=1
    return coefs[0], coefs[1]

def RE_FormFactor(magQ,ion):
    """This uses the dipole approximation.
    Note that Q must be a scalar in inverseA"""
    s = magQ/(4.*np.pi)
    
    coefs0, coefs2 = importRE_FF(ion)
    
    j0 = coefs0[0]*np.exp(-coefs0[1]*s**2) + coefs0[2]*np.exp(-coefs0[3]*s**2)  + \
        coefs0[4]*np.exp(-coefs0[5]*s**2) +coefs0[6]

    j2 = s**2*(coefs2[0]*np.exp(-coefs2[1]*s**2) + coefs2[2]*np.exp(-coefs2[3]*s**2)+\
               coefs2[4]*np.exp(-coefs2[5]*s**2) + coefs2[6])

    S, L, J = ION_NUMS_RARE_EARTH[ion]
     
    j2factor = (J*(J+1.) - S*(S+1.) + L*(L+1.))/(3.*J*(J+1.) + S*(S+1.) - L*(L+1.))
    return (j0 + j2*j2factor)**2