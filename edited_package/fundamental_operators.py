"""
Quantum mechanics operators for angular momentum calculations.

This module provides classes for angular momentum algebra in quantum mechanics,
supporting both J-basis (rare earth ions) and LS-coupling (transition metals).

Classes:
    Ket: Quantum state vector with angular momentum operations
    Operator: Matrix representation of J operators (Jx, Jy, Jz, J±)
    LSOperator: Matrix operators for LS-coupling scheme with separate L and S

Features:
    - Angular momentum ladder operators (raising/lowering)
    - Rotation operations using Euler angles and Wigner D-matrices
    - Operator arithmetic (addition, multiplication, powers)
    - Magnetization and susceptibility calculations for LS systems
    - Support for half-integer and integer angular momenta

Typical usage:
    # Create operators for j=3.5 (e.g., Yb3+)
    >>> Jz = Operator.Jz(3.5)
    >>> Jx = Operator.Jx(3.5)
    
    # Build Hamiltonian and diagonalize
    >>> H = 0.5*Jz**2 - 0.1*Jx**2
    >>> eigenvalues, eigenvectors = np.linalg.eigh(H.O)
    
    # For transition metals with LS coupling
    >>> Lz = LSOperator.Lz(L=2, S=1)
    >>> Sz = LSOperator.Sz(L=2, S=1)

References:
    - Sakurai & Napolitano, "Modern Quantum Mechanics"
    - Wigner rotation formulas (eq. 3.9.33)

Note:
    All energies in meV, magnetic fields in Tesla, temperatures in Kelvin.
"""

import numpy as np
import scipy.linalg as LA

# used inside the package
# used outside the package

class Ket():
    """
    Represents a quantum state (ket) in angular momentum basis.
    
    Provides angular momentum operators (J) and rotation operations on quantum states.
    States are represented as superpositions of |j,m⟩ basis states.
    
    Attributes:
        ket: State vector as numpy array
        j: Total angular momentum quantum number
        m: Array of magnetic quantum numbers from -j to +j
    """
    
    def __init__(self, array): 
        """
        Initialize a quantum state from coefficient array.
        
        Args:
            array: Coefficients for basis states. Length determines j via (2j+1).
                  Example: [0,0,0,1,1] represents |j=2,m=-1⟩ + |j=2,m=0⟩
        """
        self.ket = np.array(array)
        self.j = len(array)/2.0 - 0.5
        self.m = np.arange(-self.j,self.j+1,1)

    def Jz(self):
        """Apply Jz operator: returns eigenvalue m times the state."""
        return Ket( self.ket * self.m )

    def Jplus(self):
        """
        Apply J+ raising operator.
        
        J+ |j,m⟩ = √[(j-m)(j+m+1)] |j,m+1⟩
        """
        newvals = np.sqrt((self.j-self.m)*(self.j+self.m+1)) * self.ket
        return Ket( np.roll(newvals,1) )

    def Jminus(self):
        """
        Apply J- lowering operator.
        
        J- |j,m⟩ = √[(j+m)(j-m+1)] |j,m-1⟩
        """
        newvals = np.sqrt((self.j+self.m)*(self.j-self.m+1)) * self.ket
        return Ket( np.roll(newvals,-1) )

    def Jx(self):
        """
        Apply Jx operator.
        
        Jx = (J+ + J-)/2
        """
        return Ket(0.5*(self.Jplus().ket + self.Jminus().ket) )

    def Jy(self):
        """
        Apply Jy operator.
        
        Jy = -i(J+ - J-)/2
        """
        return Ket(-1j*0.5*(self.Jplus().ket - self.Jminus().ket) )

    def R(self, alpha, beta, gamma):
        """
        Rotate state using Euler angles (ZYZ convention).
        
        Args:
            alpha: First rotation about z-axis (radians)
            beta: Rotation about y-axis (radians)
            gamma: Second rotation about z-axis (radians)
            
        Returns:
            Rotated ket
        """
        return self._Rz(alpha)._Ry(beta)._Rz(gamma)

    def _Rz(self,theta):
        """
        Rotate about z-axis.
        
        Rz(θ) |j,m⟩ = exp(-imθ) |j,m⟩
        """
        newvals = np.zeros(len(self.ket), dtype=complex)
        for i in range(len(self.ket)):
           newvals[i] = self.ket[i]* np.exp(-1j*self.m[i] * theta)
        return Ket(newvals)

    def _Ry(self,beta):
        """
        Rotate about y-axis using Wigner's formula.
        
        Uses Wigner D-matrix elements for rotation.
        """
        newvals = np.zeros(len(self.ket), dtype=complex)
        for i in range(len(self.ket)):
            mm = self.m[i]
            for j in range(len(self.ket)):
                mmp = self.m[j]
                newvals[j] += self.ket[i]* self._WignersFormula(mm,mmp,beta)
        return Ket(newvals)

    def _WignersFormula(self,m,mp,beta):
        """
        Calculate Wigner d-matrix element d^j_{m,m'}(β).
        
        Args:
            m: Initial magnetic quantum number
            mp: Final magnetic quantum number
            beta: Rotation angle (radians)
            
        Returns:
            Wigner d-matrix element
            
        Reference:
            Sakurai/Napolitano eq. 3.9.33
            Cross-checked with Mathematica's WignerD function
        """
        # determine the limit of the sum over k
        kmin = np.maximum(0, m-mp)
        kmax = np.minimum(self.j+m, self.j-mp)

        d = 0
        for k in np.arange(kmin,kmax+1):
            d += (-1)**(k-m+mp) * np.sqrt(np.math.factorial(self.j+m) * np.math.factorial(self.j-m) *\
                                np.math.factorial(self.j+mp) * np.math.factorial(self.j-mp))/\
                (np.math.factorial(self.j+m-k) * np.math.factorial(k) * np.math.factorial(self.j-k-mp)*\
                 np.math.factorial(k-m+mp))*\
                np.cos(beta/2)**(2*self.j -2*k+m-mp) * np.sin(beta/2)**(2*k-m+mp)
        return d


    def __mul__(self,other):
        """
        Multiplication: inner product with ket, or matrix multiplication.
        
        Args:
            other: Another Ket (returns inner product) or matrix (returns matrix*ket)
        """
        if isinstance(other, Ket):
            # Compute inner product
            return np.dot(np.conjugate(self.ket), other.ket)
        else:
            return Ket( np.dot(self.ket, other))

    def __add__(self,other):
        """Add two kets: superposition of states."""
        if isinstance(other, Ket):
            return Ket(self.ket + other.ket)
        else:
            print("other is not a ket")


class Operator():
    """
    Angular momentum operator in matrix representation.
    
    Provides J operators (Jx, Jy, Jz, J±) as matrices in the |j,m⟩ basis.
    Supports operator arithmetic (addition, multiplication, powers).
    
    Attributes:
        O: Operator matrix (2j+1 × 2j+1)
        j: Total angular momentum quantum number
        m: Array of m values from -j to +j
    """
    
    def __init__(self, J):
        """
        Initialize operator for given j.
        """
        self.O = np.zeros((int(2*J+1), int(2*J+1)))
        self.m = np.arange(-J,J+1,1)
        self.j = J

    @staticmethod
    def Jz(J):
        """
        Create Jz operator matrix.
        
        Jz is diagonal with eigenvalues m.
        """
        obj = Operator(J)
        for i in range(len(obj.O)):
            for k in range(len(obj.O)):
                if i == k:
                    obj.O[i,k] = (obj.m[k])
        return obj

    @staticmethod
    def Jplus(J):
        """
        Create J+ raising operator matrix.
        
        J+ has off-diagonal elements √[(j-m)(j+m+1)].
        """
        obj = Operator(J)
        for i in range(len(obj.O)):
            for k in range(len(obj.O)):
                if k+1 == i:
                    obj.O[i,k] = np.sqrt((obj.j-obj.m[k])*(obj.j+obj.m[k]+1))
        return obj

    @staticmethod
    def Jminus(J):
        """
        Create J- lowering operator matrix.
        
        J- has off-diagonal elements √[(j+m)(j-m+1)].
        """
        obj = Operator(J)
        for i in range(len(obj.O)):
            for k in range(len(obj.O)):
                if k-1 == i:
                    obj.O[i,k] = np.sqrt((obj.j+obj.m[k])*(obj.j-obj.m[k]+1))
        return obj

    @staticmethod
    def Jx(J):
        """
        Create Jx operator: Jx = (J+ + J-)/2.
        """
        objp = Operator.Jplus(J)
        objm = Operator.Jminus(J)
        return 0.5*objp + 0.5*objm

    @staticmethod
    def Jy(J):
        """
        Create Jy operator: Jy = -i(J+ - J-)/2.
        """
        objp = Operator.Jplus(J)
        objm = Operator.Jminus(J)
        return -0.5j*objp + 0.5j*objm

    def __add__(self,other):
        """Add two operators or add scalar to diagonal."""
        newobj = Operator(self.j)
        if isinstance(other, Operator):
           newobj.O = self.O + other.O
        else:
           newobj.O = self.O + other*np.identity(int(2*self.j+1))
        return newobj

    def __radd__(self,other):
        """Right addition (scalar + operator)."""
        newobj = Operator(self.j)
        if isinstance(other, Operator):
            newobj.O = self.O + other.O
        else:
            newobj.O = self.O + other*np.identity(int(2*self.j+1))
        return newobj

    def __sub__(self,other):
        """Subtract operators or subtract scalar from diagonal."""
        newobj = Operator(self.j)
        if isinstance(other, Operator):
            newobj.O = self.O - other.O
        else:
            newobj.O = self.O - other*np.identity(int(2*self.j+1))
        return newobj

    def __mul__(self,other):
        """
        Multiply operator by scalar or compose two operators.
        
        Scalar multiplication or matrix multiplication O1 * O2.
        """
        newobj = Operator(self.j)
        if (isinstance(other, int) or isinstance(other, float) or isinstance(other, complex)):
           newobj.O = other * self.O
        else:
           newobj.O = np.dot(self.O, other.O)
        return newobj

    def __rmul__(self,other):
        """Right multiplication (scalar * operator or operator composition)."""
        newobj = Operator(self.j)
        if (isinstance(other, int) or isinstance(other, float)  or isinstance(other, complex)):
           newobj.O = other * self.O
        else:
           newobj.O = np.dot(other.O, self.O)
        return newobj

    def __pow__(self, power):
        """
        Raise operator to integer power.
        
        Computes O^n by repeated matrix multiplication.
        """
        newobj = Operator(self.j)
        newobj.O = self.O
        for i in range(power-1):
            newobj.O = np.dot(newobj.O,self.O)
        return newobj

    def __neg__(self):
        """Negate operator: returns -O."""
        newobj = Operator(self.j)
        newobj.O = -self.O
        return newobj

    def __repr__(self):
        """String representation shows the matrix."""
        return repr(self.O)


class LSOperator():
    """
    Operator for LS coupling (intermediate coupling scheme).
    
    Handles transition metal ions where L and S are separate good quantum numbers.
    Basis states are |L,m_L⟩ ⊗ |S,m_S⟩ with dimension (2L+1)(2S+1).
    
    Provides separate L and S operators plus arithmetic operations.
    
    Attributes:
        O: Operator matrix
        L: Orbital angular momentum quantum number
        S: Spin angular momentum quantum number
        Lm: Array of m_L values for each basis state
        Sm: Array of m_S values for each basis state
    """
    
    def __init__(self, L, S):
        """Initialize LS operator."""
        self.O = np.zeros((int((2*L+1)*(2*S+1)), int((2*L+1)*(2*S+1)) ))
        self.L = L
        self.S = S
        lm = np.arange(-L,L+1,1)
        sm = np.arange(-S,S+1,1)
        self.Lm = np.repeat(lm, len(sm))
        self.Sm = np.tile(sm, len(lm))
    
    @staticmethod
    def Lz(L, S):
        """
        Create Lz operator (acts on orbital part only).
        
        Diagonal with m_L eigenvalues.
        """
        obj = LSOperator(L, S)
        for i in range(len(obj.O)):
            for k in range(len(obj.O)):
                if i == k:
                    obj.O[i,k] = (obj.Lm[k])
        return obj

    @staticmethod
    def Lplus(L, S):
        """
        Create L+ raising operator.
        
        Raises m_L by 1, requires m_S unchanged.
        """
        obj = LSOperator(L, S)
        for i, lm1 in enumerate(obj.Lm):
            for k, lm2 in enumerate(obj.Lm):
                if (lm1 - lm2 == 1) and (obj.Sm[i] == obj.Sm[k] ):
                    obj.O[i,k] = np.sqrt((obj.L-obj.Lm[k])*(obj.L+obj.Lm[k]+1))
        return obj

    @staticmethod
    def Lminus(L, S):
        """
        Create L- lowering operator.
        
        Lowers m_L by 1, requires m_S unchanged.
        """
        obj = LSOperator(L, S)
        for i, lm1 in enumerate(obj.Lm):
            for k, lm2 in enumerate(obj.Lm):
                if (lm2 - lm1 == 1)  and (obj.Sm[i] == obj.Sm[k]):
                    obj.O[i,k] = np.sqrt((obj.L+obj.Lm[k])*(obj.L-obj.Lm[k]+1))
        return obj

    @staticmethod
    def Lx(L, S):
        """Create Lx operator: Lx = (L+ + L-)/2."""
        objp = LSOperator.Lplus(L, S)
        objm = LSOperator.Lminus(L, S)
        return 0.5*objp + 0.5*objm

    @staticmethod
    def Ly(L, S):
        """Create Ly operator: Ly = -i(L+ - L-)/2."""
        objp = LSOperator.Lplus(L, S)
        objm = LSOperator.Lminus(L, S)
        return -0.5j*objp + 0.5j*objm
    
    ##################################
    # Spin operators

    @staticmethod
    def Sz(L, S):
        """
        Create Sz operator (acts on spin part only).
        
        Diagonal with m_S eigenvalues.
        """
        obj = LSOperator(L, S)
        for i in range(len(obj.O)):
            for k in range(len(obj.O)):
                if i == k:
                    obj.O[i,k] = (obj.Sm[k])
        return obj

    @staticmethod
    def Splus(L, S):
        """
        Create S+ raising operator.
        
        Raises m_S by 1, requires m_L unchanged.
        """
        obj = LSOperator(L, S)
        for i, sm1 in enumerate(obj.Sm):
            for k, sm2 in enumerate(obj.Sm):
                if (sm1 - sm2 == 1) and (obj.Lm[i] == obj.Lm[k]):
                    obj.O[i,k] = np.sqrt((obj.S-obj.Sm[k])*(obj.S+obj.Sm[k]+1))
        return obj

    @staticmethod
    def Sminus(L, S):
        """
        Create S- lowering operator.
        
        Lowers m_S by 1, requires m_L unchanged.
        """
        obj = LSOperator(L, S)
        for i, sm1 in enumerate(obj.Sm):
            for k, sm2 in enumerate(obj.Sm):
                if (sm2 - sm1 == 1) and (obj.Lm[i] == obj.Lm[k]):
                    obj.O[i,k] = np.sqrt((obj.S+obj.Sm[k])*(obj.S-obj.Sm[k]+1))
        return obj

    @staticmethod
    def Sx(L, S):
        """Create Sx operator: Sx = (S+ + S-)/2."""
        objp = LSOperator.Splus(L, S)
        objm = LSOperator.Sminus(L, S)
        return 0.5*objp + 0.5*objm

    @staticmethod
    def Sy(L, S):
        """Create Sy operator: Sy = -i(S+ - S-)/2."""
        objp = LSOperator.Splus(L, S)
        objm = LSOperator.Sminus(L, S)
        return -0.5j*objp + 0.5j*objm


    def __add__(self,other):
        """Add operators or add scalar to diagonal."""
        newobj = LSOperator(self.L, self.S)
        try:
            newobj.O = np.add(self.O, other.O)
        except AttributeError:
            newobj.O = self.O + other*np.identity(len(self.O))
        return newobj

    def __radd__(self,other):
        """Right addition."""
        newobj = LSOperator(self.L, self.S)
        try:
            newobj.O = np.add(other.O, self.O)
        except AttributeError:
            newobj.O = self.O + other*np.identity(len(self.O))
        return newobj

    def __sub__(self,other):
        """Subtract operators or subtract scalar from diagonal."""
        newobj = LSOperator(self.L, self.S)
        try:
            newobj.O = self.O - other.O
        except AttributeError:
            newobj.O = self.O - other*np.identity(len(self.O))
        return newobj

    def __mul__(self,other):
        """Multiply by scalar or compose operators."""
        newobj = LSOperator(self.L, self.S)
        try:
            newobj.O = np.dot(self.O, other.O)
        except AttributeError:
            newobj.O = other * self.O
        return newobj

    def __rmul__(self,other):
        """Right multiplication."""
        newobj = LSOperator(self.L, self.S)
        try:
            newobj.O = np.dot(other.O, self.O)
        except AttributeError:
            newobj.O = other * self.O
        return newobj

    def __pow__(self, power):
        """Raise operator to integer power."""
        newobj = LSOperator(self.L, self.S)
        newobj.O = self.O
        for i in range(power-1):
            newobj.O = np.dot(newobj.O,self.O)
        return newobj

    def __neg__(self):
        """Negate operator."""
        newobj = LSOperator(self.L, self.S)
        newobj.O = -self.O
        return newobj

    def __repr__(self):
        """String representation shows matrix."""
        return repr(self.O)

    def magnetization(self, Temp, Field):
        """
        Calculate magnetization as function of temperature and field.
        
        Computes thermal expectation value ⟨J⟩ = Tr[ρ J] where ρ is Boltzmann
        distribution and J = L + g₀S (total angular momentum with g-factor).
        
        Args:
            Temp: Temperature(s) in Kelvin (scalar or array)
            Field: Magnetic field vector [Bx, By, Bz] in Tesla
            
        Returns:
            Magnetization vector [Mx, My, Mz] in Bohr magnetons
            If Temp is array, returns shape (3, len(Temp))
            
        Raises:
            TypeError: If Field is not 3-component vector
        """
        if len(Field) != 3: 
            raise TypeError("Field needs to be 3-component vector")

        # A) Define magnetic Hamiltonian
        Lx = LSOperator.Lx(self.L, self.S)
        Ly = LSOperator.Ly(self.L, self.S)
        Lz = LSOperator.Lz(self.L, self.S)
        Sx = LSOperator.Sx(self.L, self.S)
        Sy = LSOperator.Sy(self.L, self.S)
        Sz = LSOperator.Sz(self.L, self.S)

        g0 = 2.002319
        Jx = Lx + g0*Sx
        Jy = Ly + g0*Sy
        Jz = Lz + g0*Sz

        muB = 5.7883818012e-2  # meV/T
        JdotB = muB*((Field[0]*Lx + Field[1]*Ly + Field[2]*Lz) +\
                        (Field[0]*Sx + Field[1]*Sy + Field[2]*Sz))

        # B) Diagonalize full Hamiltonian
        FieldHam = self.O + JdotB.O
        diagonalH = LA.eigh(FieldHam)

        minE = np.amin(diagonalH[0])
        evals = diagonalH[0] - minE
        evecs = diagonalH[1].T
        # These ARE actual eigenvalues.

        # C) Compute expectation value along field
        JexpVals = np.zeros((len(evals),3))
        for i, ev in enumerate(evecs):
            JexpVals[i] =[np.real(np.dot(np.conjugate(ev), np.dot( Jx.O ,ev))),
                          np.real(np.dot(np.conjugate(ev), np.dot( Jy.O ,ev))),
                          np.real(np.dot(np.conjugate(ev), np.dot( Jz.O ,ev)))]
        k_B = 8.6173303e-2  # meV/K

        if (isinstance(Temp, int) or isinstance(Temp, float)):
            Zz = np.sum(np.exp(-evals/(k_B*Temp)))
            JexpVal = np.dot(np.exp(-evals/(k_B*Temp)),JexpVals)/Zz
            return np.real(JexpVal)
        else:
            expvals, temps = np.meshgrid(evals, Temp)
            ZZ = np.sum(np.exp(-expvals/temps/k_B), axis=1)
            JexpValList = np.repeat(JexpVals.reshape((1,)+JexpVals.shape), len(Temp), axis=0)
            JexpValList = np.sum(np.exp(-expvals/temps/k_B)*\
                                np.transpose(JexpValList, axes=[2,0,1]), axis=2) / ZZ
            return np.nan_to_num(JexpValList.T)


    def susceptibility(self, Temps, Field, deltaField):
        """
        Calculate magnetic susceptibility χ = dM/dH numerically.
        
        Uses 4-point finite difference formula for numerical derivative.
        For scalar Field, computes powder average over x, y, z directions.
        
        Args:
            Temps: Temperature array in Kelvin
            Field: Either scalar (powder average) or 3-vector (single direction)
            deltaField: Step size for numerical derivative (Tesla)
            
        Returns:
            Susceptibility in emu/mol or SI units (depending on magnetization units)
            For powder average: average of χ_xx, χ_yy, χ_zz
            For vector field: χ along field direction
            
        Raises:
            TypeError: If deltaField is not scalar
        """
        if not isinstance(deltaField, float):
            raise TypeError("Deltafield needs to be a scalar")

        if isinstance(Field, float):
            # Assume we are computing a powder average
            VecField = Field * np.array([1,0,0])
            Delta = deltaField*np.array(VecField)/Field
            Mplus1 = self.magnetization(Temps, VecField + Delta)
            Mminus1= self.magnetization(Temps, VecField - Delta)
            Mplus2 = self.magnetization(Temps, VecField + 2*Delta)
            Mminus2= self.magnetization(Temps, VecField - 2*Delta)

            dMdH_x = (8*(Mplus1 - Mminus1) - (Mplus2 - Mminus2))/(12*deltaField)

            VecField = Field * np.array([0,1,0])
            Delta = deltaField*np.array(VecField)/Field
            Mplus1 = self.magnetization(Temps, VecField + Delta)
            Mminus1= self.magnetization(Temps, VecField - Delta)
            Mplus2 = self.magnetization(Temps, VecField + 2*Delta)
            Mminus2= self.magnetization(Temps, VecField - 2*Delta)

            dMdH_y = (8*(Mplus1 - Mminus1) - (Mplus2 - Mminus2))/(12*deltaField)

            VecField = Field * np.array([0,0,1])
            Delta = deltaField*np.array(VecField)/Field
            Mplus1 = self.magnetization(Temps, VecField + Delta)
            Mminus1= self.magnetization(Temps, VecField - Delta)
            Mplus2 = self.magnetization(Temps, VecField + 2*Delta)
            Mminus2= self.magnetization(Temps, VecField - 2*Delta)

            dMdH_z = (8*(Mplus1 - Mminus1) - (Mplus2 - Mminus2))/(12*deltaField)

            return (dMdH_x[:,0]+dMdH_y[:,1]+dMdH_z[:,2])/3.

        elif len(Field) == 3:
            Delta = deltaField*np.array(Field)/np.linalg.norm(Field)
            Mplus1 = self.magnetization(Temps, Field + Delta)
            Mminus1= self.magnetization(Temps, Field - Delta)
            Mplus2 = self.magnetization(Temps, Field + 2*Delta)
            Mminus2= self.magnetization(Temps, Field - 2*Delta)

            dMdH = (8*(Mplus1 - Mminus1) - (Mplus2 - Mminus2))/(12*deltaField)

            return dMdH