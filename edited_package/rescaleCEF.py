from .constants import theta, calculate_radial_integral_RE 

# used outside the package

def rescaleCEF(ion1, ion2, B, n):
    '''Uses the point charge model to scale a CEF parameter (B)
    from one ion to another. Assumes precisely the same ligand environment
    with a different magnetic ion thrown in.
    
    Args:
        ion1: string, first ion of known CEF parameters ('Ce3+', 'Pr3+', 'Nd3+', 'Pm3+',
                                                         'Sm3+', 'Eu3+', 'Gd3+', 'Tb3+',
                                                         'Dy3+', 'Ho3+', 'Er3+', 'Tm3+', or 'Yb3+')
        ion2: string, second ion to calculate the CEF parameters for.
        B: float, a single CEF parameter for ion1.
        n: int, Operator O degree.

    Returns:
        float or int: scaled CEF parameter from Ion1 to Ion2.

    '''
    scalefact = (calculate_radial_integral_RE(ion2,n)*theta(ion2,n))/(calculate_radial_integral_RE(ion1,n)*theta(ion1,n))
    return B*scalefact