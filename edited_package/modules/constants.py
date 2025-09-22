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
JionTM = {
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
Jion = {
    'Ce3+': [0.5, 3., 2.5], 'Pr3+': [1., 5., 4.], 'Nd3+': [1.5, 6., 4.5],
    'Pm3+': [2., 6., 4.], 'Sm3+': [2.5, 5, 2.5], 'Eu3+': [3, 3, 0],
    'Gd3+': [7/2, 0, 7/2],
    'Tb3+': [3., 3., 6.], 'Dy3+': [2.5, 5., 7.5], 'Ho3+': [2., 6., 8.],
    'Er3+': [1.5, 6., 7.5], 'Tm3+': [1., 5., 6.], 'Yb3+': [0.5, 3., 3.5],
    'U4+': [1., 5., 4.], 'U3+': [1.5, 6., 4.5],
}

#######################
# Other Constants
#######################

ahc = 1.43996e4  #Constant to get the energy in units of meV = alpha*hbar*c
a0 = 0.52917721067    #Bohr radius in \AA
muB = 5.7883818012e-2  # meV/T
k_B = 8.6173303e-2  # meV/K