"""
Stevens operator equivalent matrices for crystal field calculations.

Provides Stevens operator O_n^m matrices constructed from angular momentum operators.
These are the basis operators for expanding crystal field Hamiltonians in rare earth
and transition metal systems.

Functions:
    StevensOp(J, n, m): Generate Stevens operator O_n^m for J-basis (rare earths)
    LS_StevensOp(L, S, n, m): Generate Stevens operator for LS-coupling basis (transition metals)

Stevens operators are polynomial combinations of Jz, J+, J- that transform as
spherical harmonics under rotation. They satisfy:
    - Hermiticity
    - Time-reversal symmetry: O_n^(-m) = (-1)^m (O_n^m)†
    - Tesseral form (real linear combinations)

Supported ranks: n = 0, 1, 2, 3, 4, 6 (crystal field requires n ≤ 6)
Orders: m = -n to +n

Example:
    >>> O20 = StevensOp(J=3.5, n=2, m=0)  # Quadrupole operator for Yb3+
    >>> O44 = StevensOp(J=7.5, n=4, m=4)  # Hexadecapole for Dy3+

References:
    - M.T. Hutchings, Solid State Physics 16, 227 (1964)
    - K.W.H. Stevens, Proc. Phys. Soc. A 65, 209 (1952)
"""

import numpy as np
from .fundamental_operators import Operator

# used inside the package

def StevensOp(J,n,m):
    """
    Generate stevens operator for a given total angular momentum
    and a given n and m state
    """
    Jz = Operator.Jz(J=J)
    Jp = Operator.Jplus(J=J)
    Jm = Operator.Jminus(J=J)
    X = J*(J+1.)

    if [n,m] == [0,0]:
        return np.zeros((int(2*J+1), int(2*J+1)))
    elif [n,m] == [1,0]:
        matrix = Jz
    elif [n,m] == [1,1]:
        matrix = 0.5 *(Jp + Jm)
    elif [n,m] == [1,-1]:
        matrix = -0.5j *(Jp - Jm)

    elif [n,m] == [2,2]:
        matrix = 0.5 *(Jp**2 + Jm**2)
    elif [n,m] == [2,1]:
        matrix = 0.25*(Jz*(Jp + Jm) + (Jp + Jm)*Jz)
    elif [n,m] == [2,0]:
        matrix = 3*Jz**2 - X
    elif [n,m] == [2,-1]:
        matrix = -0.25j*(Jz*(Jp - Jm) + (Jp - Jm)*Jz)
    elif [n,m] == [2,-2]:
        matrix = -0.5j *(Jp**2 - Jm**2)

    elif [n,m] == [3,3]:
        matrix = 0.5 *(Jp**3 + Jm**3)
    elif [n,m] == [3,2]:
        matrix = 0.25 *((Jp**2 + Jm**2)*Jz + Jz*(Jp**2 + Jm**2))
    elif [n,m] == [3,1]:
        matrix = 0.25*((Jp + Jm)*(5*Jz**2 - X - 0.5) + (5*Jz**2 - X - 0.5)*(Jp + Jm))
    elif [n,m] == [3,0]:
        matrix = 5*Jz**3 - (3*X-1)*Jz
    elif [n,m] == [3,-1]:
        matrix = -0.25j*((Jp - Jm)*(5*Jz**2 - X - 0.5) + (5*Jz**2 - X - 0.5)*(Jp - Jm))
    elif [n,m] == [3,-2]:
        matrix = -0.25j*(Jz*(Jp**2 - Jm**2) + (Jp**2 - Jm**2)*Jz)
    elif [n,m] == [3,-3]:
        matrix = -0.5j *(Jp**3 - Jm**3)

    elif [n,m] == [4,4]:
        matrix = 0.5 *(Jp**4 + Jm**4)
    elif [n,m] == [4,3]:
        matrix = 0.25 *((Jp**3 + Jm**3)*Jz + Jz*(Jp**3 + Jm**3))
    elif [n,m] == [4,2]:
        matrix = 0.25 *((Jp**2 + Jm**2)*(7*Jz**2 -X -5) + (7*Jz**2 -X -5)*(Jp**2 + Jm**2))
    elif [n,m] == [4,1]:
        matrix = 0.25 *((Jp + Jm)*(7*Jz**3 -(3*X+1)*Jz) + (7*Jz**3 -(3*X+1)*Jz)*(Jp + Jm))
    elif [n,m] == [4,0]:
        matrix = 35*Jz**4 - (30*X -25)*Jz**2 + 3*X**2 - 6*X
    elif [n,m] == [4,-4]:
        matrix = -0.5j *(Jp**4 - Jm**4)
    elif [n,m] == [4,-3]:
        matrix = -0.25j *((Jp**3 - Jm**3)*Jz + Jz*(Jp**3 - Jm**3))
    elif [n,m] == [4,-2]:
        matrix = -0.25j *((Jp**2 - Jm**2)*(7*Jz**2 -X -5) + (7*Jz**2 -X -5)*(Jp**2 - Jm**2))
    elif [n,m] == [4,-1]:
        matrix = -0.25j *((Jp - Jm)*(7*Jz**3 -(3*X+1)*Jz) + (7*Jz**3 -(3*X+1)*Jz)*(Jp - Jm))

    elif [n,m] == [6,6]:
        matrix = 0.5 *(Jp**6 + Jm**6)
    elif [n,m] == [6,5]:
        matrix = 0.25*((Jp**5 + Jm**5)*Jz + Jz*(Jp**5 + Jm**5))
    elif [n,m] == [6,4]:
        matrix = 0.25*((Jp**4 + Jm**4)*(11*Jz**2 -X -38) + (11*Jz**2 -X -38)*(Jp**4 + Jm**4))
    elif [n,m] == [6,3]:
        matrix = 0.25*((Jp**3 + Jm**3)*(11*Jz**3 -(3*X+59)*Jz) + (11*Jz**3 -(3*X+59)*Jz)*(Jp**3 + Jm**3))
    elif [n,m] == [6,2]:
        matrix = 0.25*((Jp**2 + Jm**2)*(33*Jz**4 -(18*X+123)*Jz**2 +X**2 +10*X +102) +\
                    (33*Jz**4 -(18*X+123)*Jz**2 +X**2 +10*X +102)*(Jp**2 + Jm**2))
    elif [n,m] == [6,1]:
        matrix = 0.25*((Jp +Jm)*(33*Jz**5 -(30*X-15)*Jz**3 +(5*X**2 -10*X +12)*Jz) +\
                    (33*Jz**5 -(30*X-15)*Jz**3 +(5*X**2 -10*X +12)*Jz)*(Jp+ Jm))
    elif [n,m] == [6,0]:
        matrix = 231*Jz**6 - (315*X-735)*Jz**4 + (105*X**2 -525*X +294)*Jz**2 -\
                 5*X**3 + 40*X**2 - 60*X
    elif [n,m] == [6,-6]:
        matrix = -0.5j *(Jp**6 - Jm**6)
    elif [n,m] == [6,-5]:
        matrix = -0.25j*((Jp**5 - Jm**5)*Jz + Jz*(Jp**5 - Jm**5))
    elif [n,m] == [6,-4]:
        matrix = -0.25j*((Jp**4 - Jm**4)*(11*Jz**2 -X -38) + (11*Jz**2 -X -38)*(Jp**4 - Jm**4))
    elif [n,m] == [6,-3]:
        matrix = -0.25j*((Jp**3 - Jm**3)*(11*Jz**3 -(3*X+59)*Jz) + (11*Jz**3 -(3*X+59)*Jz)*(Jp**3 - Jm**3))
    elif [n,m] == [6,-2]:
        matrix = -0.25j*((Jp**2 - Jm**2)*(33*Jz**4 -(18*X+123)*Jz**2 +X**2 +10*X +102) +\
                    (33*Jz**4 -(18*X+123)*Jz**2 +X**2 +10*X +102)*(Jp**2 - Jm**2))
    elif [n,m] == [6,-1]:
        matrix = -0.25j*((Jp - Jm)*(33*Jz**5 -(30*X-15)*Jz**3 +(5*X**2 -10*X +12)*Jz) +\
                 (33*Jz**5 -(30*X-15)*Jz**3 +(5*X**2 -10*X +12)*Jz)*(Jp - Jm))

    return matrix.O


def LS_StevensOp(L,S,n,m):
    """generate stevens operator for a given total angular momentum
    and a given n and m state, but in the LS basis"""
    lmatrix = StevensOp(L,n,m)

    fullmatrix = np.hstack(np.hstack(np.multiply.outer(lmatrix, np.identity(int(2*S+1)))))
    return fullmatrix
