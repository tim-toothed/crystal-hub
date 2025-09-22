from constants import Jion

#######################
# LandeGFactor Function
#######################

def LandeGFactor(ion):
    """
    Calculate the Lande g-factor for a given ion. 
    Lande g-factor describes the magnetic moment of an atom or ion.
    
    Args:
        ion (str): Ion name (e.g., 'Ce3+', 'Pr3+')
        
    Returns:
        float: Lande g-factor value
    """
    s, l, j = Jion[ion]
    return 1.5 + (s*(s+1.) - l*(l+1.))/(2.*j*(j+1.))