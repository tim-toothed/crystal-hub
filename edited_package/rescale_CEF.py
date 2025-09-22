from pcf_lib.PointChargeConstants import theta, RadialIntegral 

def rescaleCEF(ion1, ion2, B, n):
    '''Uses the point charge model to scale a CEF parameter (B)
    from one ion to another. Assumes precisely the same ligand environment
    with a different magnetic ion thrown in.'''
    scalefact = (RadialIntegral(ion2,n)*theta(ion2,n))/(RadialIntegral(ion1,n)*theta(ion1,n))
    return B*scalefact