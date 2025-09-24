def IsHalfFilled(ion):
    '''determine if given ion has a half-filled shell or not.'''

    HalfList = ['Mn2+',
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
    notHalfList = ['Ti2+','Ti3+',
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

    if ion in HalfList:
        return True
    elif ion in notHalfList:
        return False
    else:
        raise ValueError('{} is not a known ion for PyCrystalField.'.format(ion))