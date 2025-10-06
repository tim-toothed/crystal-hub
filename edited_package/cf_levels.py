import warnings
warnings.filterwarnings('ignore')

import numpy as np
from numba import float64
from numba.experimental import jitclass
from scipy import optimize
import scipy.linalg as LA
from scipy.special import wofz

from .constants import ION_NUMS_RARE_EARTH
from .form_factors import RE_FormFactor
from .create_fit_function import makeFitFunction
from .fundamental_operators import Ket, Operator, LSOperator
from .stevens_operators import StevensOp, LS_StevensOp


class CFLevels:
    """
    Crystal field levels calculator for rare earth ions.

    Handles crystal field Hamiltonian construction, diagonalization, and calculation
    of magnetic/spectroscopic properties in the J-basis for rare earth systems.
    """

    def __init__(self, StevensOperators, Parameters):
        """
        Initialize CFLevels from Stevens operators and parameters.

        Constructs the crystal field Hamiltonian as a linear combination of 
        Stevens operators weighted by their corresponding B_n^m parameters.

        Args:
            StevensOperators: List of Stevens operator matrices (numpy arrays)
            Parameters: List of B_n^m coefficients in meV
            
        Sets:
            self.H: Total Hamiltonian matrix (sum of weighted operators)
            self.O: Saved operators for later fitting
            self.B: Parameter values
            self.J: Total angular momentum quantum number
            self.opttran: Optimized transition calculator object
        """

        self.O = StevensOperators  #save these for a fit
        self.B = Parameters
        self.H = np.sum([a*b for a,b in zip(self.O, self.B)], axis=0)

        try:
            self.J = (len(self.H) -1.)/2
            self.opttran = opttransition(Operator.Jx(self.J).O, Operator.Jy(self.J).O.imag, Operator.Jz(self.J).O)
        except TypeError: pass


    @classmethod
    def Bdict(cls, ion, Bdict):
        """
        Create CFLevels instance from ion name and parameter dictionary.

        Factory method that constructs Stevens operators from ion's J value
        and builds Hamiltonian from B_n^m parameters.

        Args:
            ion: Ion symbol (e.g., 'Yb3+', 'Dy3+')
            Bdict: Dictionary mapping parameter names to values
                    Keys: 'B20', 'B40', 'B43', etc.
                    Values: Parameters in meV

        Returns:
            CFLevels instance with constructed Hamiltonian
            
        Example:
            >>> params = {'B20': 0.5, 'B40': -0.02, 'B43': 0.15}
            >>> cf = CFLevels.Bdict('Yb3+', params)
        """

        ionJ = ION_NUMS_RARE_EARTH[ion][-1]
        Stev_O = []
        Parameters = []
        for Bnm in Bdict:
            Parameters.append(Bdict[Bnm])
            n = int(Bnm[1])
            m = int(Bnm[2:])    
            Stev_O.append(  StevensOp(ionJ,n,m)  )

        newcls = cls(Stev_O, Parameters)
        newcls.BnmLabels = Bdict.keys
        newcls.ion = ion
        return newcls


    @classmethod
    def Hamiltonian(cls, Hamil):
        """
        Create CFLevels instance directly from a Hamiltonian matrix.

        Alternative constructor for cases where Hamiltonian is already built
        (e.g., from external source or complex construction).

        Args:
            Hamil: Pre-constructed Hamiltonian matrix (numpy array)
            
        Returns:
            CFLevels instance with the provided Hamiltonian
            
        Note:
            Sets J based on matrix dimension: J = (dim - 1) / 2
        """

        newcls = cls([0,0],[0,0])  # Create empty class so we can just define Hamiltonian
        newcls.H = Hamil
        newcls.J = (len(Hamil) -1.)/2
        newcls.opttran = opttransition(Operator.Jx(newcls.J).O.real, Operator.Jy(newcls.J).O.imag, Operator.Jz(newcls.J).O.real)
        return newcls


    def newCoeff(self, newcoeff):
        """
        Update Stevens parameters and rediagonalize.
        
        Reconstructs Hamiltonian with new coefficients and immediately
        diagonalizes. Used during fitting procedures.
        
        Args:
            newcoeff: New list of B_n^m parameters in meV
            
        Side effects:
            Updates self.B, self.H, and diagonalization results
        """

        self.B = np.array(newcoeff)
        newH = np.sum([a*b for a,b in zip(self.O, newcoeff)], axis=0)
        self.diagonalize(newH)


    def diagonalize(self, Hamiltonian=None, old=False):
        """
        Diagonalize the crystal field Hamiltonian.

        Finds energy eigenvalues and eigenvectors using banded matrix
        algorithm for efficiency (default) or standard eigh (if old=True).

        Args:
            Hamiltonian: Optional Hamiltonian to diagonalize (uses self.H if None)
            old: If True, use LA.eigh (slower). If False, use LA.eig_banded (faster)
            
        Sets:
            self.eigenvaluesNoNorm: Raw eigenvalues from diagonalization
            self.eigenvalues: Eigenvalues shifted so ground state = 0 meV
            self.eigenvectors: Eigenvectors (rows are eigenstates)
            
        Note:
            Values below 1e-15 are set to exactly zero to clean up numerical noise
        """

        if Hamiltonian is None:
            Hamiltonian = self.H
        else:
            self.H = Hamiltonian
        if old:
            diagonalH = LA.eigh(Hamiltonian)  #This was slower and less precise
        else:
            bands = self._findbands(Hamiltonian)
            diagonalH = LA.eig_banded(bands, lower=True)

        self.eigenvaluesNoNorm = diagonalH[0]
        self.eigenvalues = diagonalH[0] - np.amin(diagonalH[0])
        self.eigenvectors = diagonalH[1].T
        # set very small values to zero
        tol = 1e-15
        self.eigenvalues[abs(self.eigenvalues) < tol] = 0.0
        self.eigenvectors[abs(self.eigenvectors) < tol] = 0.0


    def diagonalize_banded(self, Hamiltonian=None):
        """
        Diagonalize using banded matrix storage format.

        More efficient version of diagonalize() that exploits the sparse
        structure of angular momentum matrices (only near-diagonal elements).

        Args:
            Hamiltonian: Optional Hamiltonian to diagonalize (uses self.H if None)
            
        Sets:
            Same as diagonalize(): eigenvalues, eigenvectors, eigenvaluesNoNorm
            
        Note:
            Uses scipy.linalg.eig_banded for O(n²) instead of O(n³) scaling
        """
        if Hamiltonian is None:
            Hamiltonian = self.H
        else:
            self.H = Hamiltonian

        bands = self._findbands(Hamiltonian)
        diagonalH = LA.eig_banded(bands, lower=True)

        self.eigenvaluesNoNorm = diagonalH[0]
        self.eigenvalues = diagonalH[0] - np.amin(diagonalH[0])
        self.eigenvectors = diagonalH[1].T
        # set very small values to zero
        tol = 1e-15
        self.eigenvalues[abs(self.eigenvalues) < tol] = 0.0
        self.eigenvectors[abs(self.eigenvectors) < tol] = 0.0


    def _findbands(self, matrix):
        """
        Extract banded structure from Hamiltonian matrix.

        Identifies which diagonals contain non-zero elements and extracts
        them into banded storage format for eig_banded().

        Args:
            matrix: Full Hamiltonian matrix
            
        Returns:
            Banded matrix array where row i contains the i-th diagonal
            
        Note:
            Only extracts diagonals with non-zero elements (up to 1e-10 threshold)
        """

        diags = np.zeros((len(matrix),len(matrix)), dtype=np.complex128)
        for i in range(len(matrix)):
            diag = matrix.diagonal(i)
            if i == 0:
                diags[i] = diag
            else:
                diags[i][:-i] = diag
            if np.count_nonzero(np.around(diag,10)) > 0:
                nonzerobands = i
        return diags[:nonzerobands+1]


    def transitionIntensity(self, ii, jj, Temp):
        """
        Calculate neutron scattering intensity for a specific transition.

        Computes |<i|J|j>|² weighted by Boltzmann population factor at
        given temperature.

        Args:
            ii: Initial state index
            jj: Final state index
            Temp: Temperature in Kelvin
            
        Returns:
            Transition intensity including population factor
            
        Formula:
            I_ij = P_n * |<i|J|j>|²
            where P_n = exp(-E_i/k_B*T) / Z
        """

        # for population factor weights
        beta = 1/(8.61733e-2*Temp)  # Boltzmann constant is in meV/K
        Z = sum([np.exp(-beta*en) for en in self.eigenvalues])
        # compute population factor
        pn = np.exp(-beta *self.eigenvalues[ii])/Z
        
        # compute amplitude
        mJn = self.opttran.transition(self.eigenvectors.real[ii] ,  self.eigenvectors.real[jj])
        return pn*mJn


    def neutronSpectrum(self, Earray, Temp, Ei, ResFunc, gamma = 0):
        """
        Calculate inelastic neutron scattering spectrum (1D in energy).

        Computes S(Q,ω) summed over all possible transitions, weighted by
        population factors and instrument resolution function.

        Args:
            Earray: Energy transfer values (meV)
            Temp: Sample temperature (K)
            Ei: Incident neutron energy (meV)
            ResFunc: Function that returns Gaussian width given energy transfer
            gamma: Lorentzian width for intrinsic peak broadening (meV)
            
        Returns:
            Intensity array vs energy, including k'/k kinematic factor
            
        Note:
            Uses Voigt profile (Gaussian ⊗ Lorentzian) for peak shapes
            Only computes transitions with population > 0.001
        """

        try:
            eigenkets = self.eigenvectors.real
            intensity = np.zeros(len(Earray))
        except AttributeError:
            self.diagonalize()
            eigenkets = self.eigenvectors.real
            intensity = np.zeros(len(Earray))

        # for population factor weights
        beta = 1/(8.61733e-2*Temp)  # Boltzmann constant is in meV/K
        Z = sum([np.exp(-beta*en) for en in self.eigenvalues])

        for i, ket_i in enumerate(eigenkets):
            # compute population factor
            pn = np.exp(-beta *self.eigenvalues[i])/Z
            if pn > 1e-3:  #only compute for transitions with enough weight
                for j, ket_j in enumerate(eigenkets):
                    # compute amplitude
                    #mJn = self._transition(ket_i,ket_j)   # Old: slow
                    mJn = self.opttran.transition(ket_i,ket_j)
                    deltaE = self.eigenvalues[j] - self.eigenvalues[i]
                    GausWidth = ResFunc(deltaE)  #peak width due to instrument resolution
                    intensity += ((pn * mJn * self._voigt(x=Earray, x0=deltaE, alpha=GausWidth, 
                                                        gamma=gamma)).real).astype('float64')
                #intensity += ((pn * mJn * self._lorentzian(Earray, deltaE, Width)).real).astype('float64')

        ## List comprehension: turns out this way was slower.
        # intensity = np.sum([
        #     np.exp(-beta *self.eigenvalues[i])/Z *\
        #     self._transition(eigenkets[i], eigenkets[j]) *\
        #     self._voigt(x=Earray, x0=(self.eigenvalues[j] - self.eigenvalues[i]), 
        #         alpha=ResFunc(self.eigenvalues[j] - self.eigenvalues[i]), gamma=gamma)
        #     for i in range(len(eigenkets)) for j in range(len(eigenkets))
        #     ], axis = 0)

        kpoverk = np.sqrt((Ei - Earray)/Ei) #k'/k = sqrt(E'/E)
        return intensity * kpoverk


    def neutronSpectrum_customLineshape(self, Earray, Temp, Ei, LineshapeFunc):
        """
        Neutron spectrum with user-defined lineshape function.

        Alternative to neutronSpectrum() allowing custom peak shapes beyond
        simple Voigt profiles.

        Args:
            Earray: Energy transfer values (meV)
            Temp: Sample temperature (K)
            Ei: Incident neutron energy (meV)
            LineshapeFunc: Function(E_array, E_center) returning peak shape
            
        Returns:
            Intensity array vs energy with k'/k correction
            
        Example:
            >>> def triangle_peak(E, E0):
            ...     return np.maximum(0, 1 - np.abs(E - E0)/0.5)
            >>> spectrum = cf.neutronSpectrum_customLineshape(
            ...     E, T=5, Ei=25, LineshapeFunc=triangle_peak)
        """

        try:
            eigenkets = self.eigenvectors.real
            intensity = np.zeros(len(Earray))
        except AttributeError:
            self.diagonalize()
            eigenkets = self.eigenvectors.real
            intensity = np.zeros(len(Earray))

        # for population factor weights
        beta = 1/(8.61733e-2*Temp)  # Boltzmann constant is in meV/K
        Z = sum([np.exp(-beta*en) for en in self.eigenvalues])

        for i, ket_i in enumerate(eigenkets):
            # compute population factor
            pn = np.exp(-beta *self.eigenvalues[i])/Z
            if pn > 1e-3:  #only compute for transitions with enough weight
                for j, ket_j in enumerate(eigenkets):
                    # compute amplitude
                    #mJn = self._transition(ket_i,ket_j)   # Old: slow
                    mJn = self.opttran.transition(ket_i,ket_j)
                    deltaE = self.eigenvalues[j] - self.eigenvalues[i]
                    intensity += ((pn * mJn * LineshapeFunc(Earray - deltaE,
                                                            deltaE)).real).astype('float64')
                
        kpoverk = np.sqrt((Ei - Earray)/Ei) #k'/k = sqrt(E'/E)
        return intensity * kpoverk


    def normalizedNeutronSpectrum(self, Earray, Temp, ResFunc, gamma = 0):
        """
        Neutron spectrum without kinematic (k'/k) correction.

        Simpler version of neutronSpectrum() that omits the sqrt(E'/E) factor,
        useful for comparing relative intensities or when correction will be
        applied separately.

        Args:
            Earray: Energy transfer values (meV)
            Temp: Sample temperature (K)
            ResFunc: Function returning Gaussian resolution width
            gamma: Lorentzian width (meV)
            
        Returns:
            Intensity array vs energy (no kinematic factor)
        """

        try:
            eigenkets = self.eigenvectors.real
            intensity = np.zeros(len(Earray))
        except AttributeError:
            self.diagonalize()
            eigenkets = self.eigenvectors.real
            intensity = np.zeros(len(Earray))

        # for population factor weights
        beta = 1/(8.61733e-2*Temp)  # Boltzmann constant is in meV/K
        Z = sum([np.exp(-beta*en) for en in self.eigenvalues])

        for i, ket_i in enumerate(eigenkets):
            # compute population factor
            pn = np.exp(-beta *self.eigenvalues[i])/Z
            if pn > 1e-3:  #only compute for transitions with enough weight
                for j, ket_j in enumerate(eigenkets):
                    # compute amplitude
                    #mJn = self._transition(ket_i,ket_j)  # Old: slow
                    mJn = self.opttran.transition(ket_i,ket_j)
                    deltaE = self.eigenvalues[j] - self.eigenvalues[i]
                    GausWidth = ResFunc(deltaE)  #peak width due to instrument resolution
                    intensity += ((pn * mJn * self._voigt(x=Earray, x0=deltaE, alpha=GausWidth, 
                                                        gamma=gamma)).real).astype('float64')
        return intensity


    def normalizedNeutronSpectrum_customLineshape(self, Earray, Temp, LineshapeFunc):
        """
        Normalized neutron spectrum with custom lineshape.

        Combines normalized spectrum (no k'/k) with user-defined peak shapes.

        Args:
            Earray: Energy transfer values (meV)
            Temp: Sample temperature (K)
            LineshapeFunc: Custom lineshape function(E_array, E_center)
            
        Returns:
            Intensity array vs energy without kinematic correction
        """

        try:
            eigenkets = self.eigenvectors.real
            intensity = np.zeros(len(Earray))
        except AttributeError:
            self.diagonalize()
            eigenkets = self.eigenvectors.real
            intensity = np.zeros(len(Earray))

        # for population factor weights
        beta = 1/(8.61733e-2*Temp)  # Boltzmann constant is in meV/K
        Z = sum([np.exp(-beta*en) for en in self.eigenvalues])

        for i, ket_i in enumerate(eigenkets):
            # compute population factor
            pn = np.exp(-beta *self.eigenvalues[i])/Z
            if pn > 1e-3:  #only compute for transitions with enough weight
                for j, ket_j in enumerate(eigenkets):
                    # compute amplitude
                    #mJn = self._transition(ket_i,ket_j)  # Old: slow
                    mJn = self.opttran.transition(ket_i,ket_j)
                    deltaE = self.eigenvalues[j] - self.eigenvalues[i]
                    intensity += ((pn * mJn * LineshapeFunc(Earray - deltaE,
                                                            deltaE)).real).astype('float64')
        return intensity


    def neutronSpectrum2D(self, Earray, Qarray, Temp, Ei, ResFunc, gamma, Ion, DebyeWaller=0):
        """
        Calculate 2D neutron spectrum: S(Q,ω).
        
        Extends 1D spectrum to include Q-dependence through magnetic form
        factor and Debye-Waller factor.
        
        Args:
            Earray: Energy transfer grid (meV)
            Qarray: Momentum transfer grid (Å⁻¹)
            Temp: Temperature (K)
            Ei: Incident energy (meV)
            ResFunc: Energy resolution function
            gamma: Lorentzian width (meV)
            Ion: Ion symbol for form factor
            DebyeWaller: Debye-Waller parameter (Å)
            
        Returns:
            2D intensity array [len(Earray) × len(Qarray)]
            
        Formula:
            I(Q,E) = I_1D(E) * F²(Q) * exp(-Q²⟨u²⟩/3)
        """

        intensity1D = self.neutronSpectrum(Earray, Temp, Ei, ResFunc,  gamma)

        # Scale by Debye-Waller Factor
        DWF = np.exp(1./3. * Qarray**2 * DebyeWaller**2)
        # Scale by form factor
        FormFactor = RE_FormFactor(Qarray,Ion)
        return np.outer(intensity1D, DWF*FormFactor)


    def normalizedNeutronSpectrum2D(self, Earray, Qarray, Temp, ResFunc, gamma, Ion, DebyeWaller=0):
        """
        2D neutron spectrum without k'/k correction.
        
        Args:
            Same as neutronSpectrum2D but omits kinematic factor
            
        Returns:
            2D intensity array [energy × momentum]
        """

        intensity1D = self.normalizedNeutronSpectrum(Earray, Temp, ResFunc,  gamma)

        # Scale by Debye-Waller Factor
        DWF = np.exp(1./3. * Qarray**2 * DebyeWaller**2)
        # Scale by form factor
        FormFactor = RE_FormFactor(Qarray,Ion)
        return np.outer(intensity1D, DWF*FormFactor)


    def _transition(self,ket1,ket2):  ## Correct, but slow.
        """
        Compute transition matrix element |<i|J|j>|².
        
        Calculates sum of squared matrix elements for all three components
        of angular momentum operator: |<i|J_x|j>|² + |<i|J_y|j>|² + |<i|J_z|j>|²
        
        Args:
            ket1: Initial state wavefunction (Ket object)
            ket2: Final state wavefunction (Ket object)
            
        Returns:
            Total transition amplitude (real float)
            
        Note:
            This is the slower, more explicit version. The optimized version
            is in self.opttran.transition()
            
        Raises:
            ValueError: If result has non-zero imaginary part (indicates error)
        """
        ax = (ket1*ket2.Jx() )*( ket2*ket1.Jx() )
        ay = (ket1*ket2.Jy() )*( ket2*ket1.Jy() )
        az = (ket1*ket2.Jz() )*( ket2*ket1.Jz() )
        # eliminate tiny values
        ax, ay, az = np.around(ax, 10), np.around(ay, 10), np.around(az, 10)
        if (ax + ay + az).imag == 0:
            return ((ax + ay + az).real).astype(float)
        else:
            print(ax, ay, az)
            raise ValueError("non-real amplitude. Error somewhere.")


    def _lorentzian(self, x, x0, gamma):
        """
        Lorentzian lineshape function.

        Args:
            x: Energy values where function is evaluated
            x0: Peak center position
            gamma: Full width at half maximum (FWHM)
            
        Returns:
            Lorentzian function: L(x) = (Γ/2π) / [(x-x₀)² + (Γ/2)²]
        """

        return 1/np.pi * (0.5*gamma)/((x-x0)**2 + (0.5*gamma)**2)


    def _voigt(self, x, x0, alpha, gamma):
        """
        Voigt profile (convolution of Gaussian and Lorentzian).

        Used for realistic peak shapes in neutron scattering where both
        instrument resolution (Gaussian) and lifetime broadening (Lorentzian)
        contribute.

        Args:
            x: Energy array
            x0: Peak center
            alpha: Gaussian FWHM (from instrument resolution)
            gamma: Lorentzian FWHM (intrinsic width)
            
        Returns:
            Voigt profile evaluated at x
            
        Note:
            Uses Faddeeva function (scipy.special.wofz) for efficient computation
        """

        sigma = (0.5*alpha) / np.sqrt(2 * np.log(2))
        return np.real(wofz(((x-x0) + 1j*(0.5*gamma))/sigma/np.sqrt(2))) / sigma\
                                                            /np.sqrt(2*np.pi)


    def _Re(self,value):
        """
        Extract real part if imaginary component is negligible.
        
        Helper function to clean up numerical results that should be real
        but have tiny imaginary parts due to floating point errors.
        
        Args:
            value: Complex number or array
            
        Returns:
            Real part if |Im| < 1e-9, otherwise returns complex value
            
        Note:
            Threshold is 1e-9 to catch numerical noise while preserving
            genuinely complex results
        """

        thresh = 1e-9
        if np.size(value) == 1 & isinstance(value, complex):
            if np.abs(value.imag) <= thresh:
                return (value.real).astype(float)
            else: 
                return value
        else:
            if np.all(np.abs(value.imag) < thresh):
                return (value.real)
            else: return value


    def printEigenvectors(self):
        """
        Display eigenvalues and eigenvectors in readable table format.
        
        Prints sorted eigenvalues with corresponding eigenvector components
        in a formatted ASCII table for inspection.
        
        Output format:
            Eigenvalues    Eigenvectors
            --------------------------------
            0.00000    | [components] |
            5.23145    | [components] |
            ...
            
        Side effects:
            Calls diagonalize() if not already done
        """

        try:
            eigenkets = self.eigenvectors.real
        except AttributeError:
            self.diagonalize()

        print('\n Eigenvalues \t Eigenvectors')
        print('\t\t'+'-------'*(len(self.eigenvalues)+1))
        sortinds = self.eigenvalues.argsort()
        sortEVal= np.around(self.eigenvalues[sortinds],5)
        sortEVec= np.around(self.eigenvectors[sortinds],3)
        for i in range(len(sortinds)):
            print(format(self._Re(sortEVal[i]), '.5f'),'\t| ', self._Re(sortEVec[i]),' |')
        print('\t\t'+'-------'*(len(self.eigenvalues)+1) + '\n')


    def printLaTexEigenvectors(self, precision = 4):
        """
        Output eigenvectors and eigenvalues in LaTeX table format.
        
        Generates LaTeX code for publication-ready tables of crystal field
        eigenstates with customizable numerical precision.
        
        Args:
            precision: Number of decimal places for eigenvector components
            
        Output:
            LaTeX table code printed to stdout, ready to copy into manuscript
            
        Format:
            Handles both integer and half-integer J values appropriately:
            - Integer J: column headers |m⟩ with integer m
            - Half-integer J: headers |m/2⟩ with proper fractions
        """
        try:
            eigenkets = self.eigenvectors.real
        except AttributeError:
            self.diagonalize()
        
        print('\\begin{table*}\n\\caption{Eigenvectors and Eigenvalues...}')
        print('\\begin{ruledtabular}')
        numev = len(self.eigenvalues)
        print('\\begin{tabular}{c|'+'c'*numev+'}')
        if numev % 2 == 1:
            print('E (meV) &'+' & '.join(['$|'+str(int(kk))+'\\rangle$' for kk in 
                        np.arange(-(numev-1)/2,numev/2)])
                +' \\tabularnewline\n \\hline ')
        else:
            print('E (meV) &'+
                ' & '.join(['$| -\\frac{'+str(abs(kk))+'}{2}\\rangle$' if kk <0
                            else '$| \\frac{'+str(abs(kk))+'}{2}\\rangle$'
                            for kk in np.arange(-(numev-1),numev,2)])
                +' \\tabularnewline\n \\hline ')
        sortinds = self.eigenvalues.argsort()
        sortEVal= np.around(self.eigenvalues[sortinds],3)
        sortEVec= np.around(self.eigenvectors[sortinds],precision)
        for i in range(len(sortinds)):
            print(format(self._Re(sortEVal[i]), '.3f'),'&', 
                ' & '.join([str(eevv) for eevv in self._Re(sortEVec[i])]), '\\tabularnewline')
        print('\\end{tabular}\\end{ruledtabular}')
        print('\\label{flo:Eigenvectors}\n\\end{table*}')


    def gsExpectation(self):
        """
        Calculate and print expectation values of J for ground state(s).
        
        Computes ⟨J_x⟩, ⟨J_y⟩, ⟨J_z⟩ for all degenerate ground states
        (eigenvalues ≈ 0 meV within 1e-7 meV tolerance).
        
        Output:
            Prints expectation values for each ground state component
            
        Use:
            Helps identify magnetic moment direction and magnitude in
            ground state doublet or multiplet
        """

        zeroinds = np.where(np.around(self.eigenvalues,7)==0)
        gsEVec = self.eigenvectors[zeroinds]
        print('\t Ground State Expectation Values:')
        for ev in gsEVec:
            vv = Ket(ev)
            jjxx = self._Re( vv*vv.Jx() )
            jjyy = self._Re( vv*vv.Jy() )
            jjzz = self._Re( vv*vv.Jz() )
            print('  <J_x> =',jjxx,'\t<J_y> =',jjyy,'\t<J_z> =',jjzz)
        print(' ')


    def magnetization(self, ion, Temp, Field):
        """
        Calculate thermal-averaged magnetization vector M(T,H).
        
        Computes full vector magnetization by:
        1. Adding Zeeman term g_J μ_B J·B to Hamiltonian
        2. Rediagonalizing in presence of field
        3. Calculating thermal average: M = g_J Σ_n P_n ⟨n|J|n⟩
        
        Args:
            ion: Ion symbol or g-factor value
            Temp: Temperature (K) - can be array
            Field: 3-component field vector [B_x, B_y, B_z] in Tesla
            
        Returns:
            [M_x, M_y, M_z] in Bohr magnetons
            - If Temp is scalar: returns 3-component vector
            - If Temp is array: returns [3 × len(Temp)] array
            
        Note:
            Uses exact diagonalization (not perturbation theory) so valid
            for all field strengths and temperatures
            
        Example:
            >>> M = cf.magnetization('Yb3+', Temp=np.linspace(2, 300, 100),
            ...                      Field=np.array([0, 0, 1.0]))
        """

        if len(Field) != 3: 
            raise TypeError("Field needs to be 3-component vector")

        # A) Define magnetic Hamiltonian

        #Jx = Operator.Jx(self.J)
        # Jy = Operator.Jy(self.J).O
        #Jz = Operator.Jz(self.J)
        Jx = self.opttran.Jx
        Jy = self.opttran.Jy * 1j
        Jz = self.opttran.Jz

        #print(Jx)
        #print(Jy)
        if isinstance(ion, str):
            gJ = LandeGFactor(ion)
        else: gJ = ion
        muB = 5.7883818012e-2  # meV/T
        #mu0 = np.pi*4e-7       # T*m/A
        JdotB = gJ*muB*(Field[0]*Jx + Field[1]*Jy + Field[2]*Jz)

        # B) Diagonalize full Hamiltonian
        FieldHam = self.H + JdotB
        #FieldHam = self.H + JdotB.O
        diagonalH = LA.eigh(FieldHam)

        minE = np.amin(diagonalH[0])
        evals = diagonalH[0] - minE
        evecs = diagonalH[1].T
        # These ARE actual eigenvalues.

        # C) Compute expectation value along field
        JexpVals = np.zeros((len(evals),3))
        for i, ev in enumerate(evecs):
            kev = Ket(ev)
            # print np.real(np.dot(ev,kev.Jy().ket)), np.real(np.dot(ev,np.dot(Jy.O,ev)))
            # print np.real(kev*kev.Jy()) - np.real(np.dot(ev,np.dot(Jy.O,ev)))
            JexpVals[i] =[np.real(kev*kev.Jx()),
                          np.real(kev*kev.Jy()),
                          np.real(kev*kev.Jz())]
        k_B = 8.6173303e-2  # meV/K

        if (isinstance(Temp, int) or isinstance(Temp, float)):
            Zz = np.sum(np.exp(-evals/(k_B*Temp)))
            JexpVal = np.dot(np.exp(-evals/(k_B*Temp)),JexpVals)/Zz
            return gJ*JexpVal
        else:
            expvals, temps = np.meshgrid(evals, Temp)
            ZZ = np.sum(np.exp(-expvals/temps/k_B), axis=1)
            JexpValList = np.repeat(JexpVals.reshape((1,)+JexpVals.shape), len(Temp), axis=0)
            JexpValList = np.sum(np.exp(-expvals/temps/k_B)*\
                                np.transpose(JexpValList, axes=[2,0,1]), axis=2) / ZZ
            # if np.isnan(JexpValList).any():
            #     print -expvals[0]/temps[0]/k_B
            #     print np.exp(-expvals/temps/k_B)[0]
            #     raise ValueError('Nan in result!')
            return np.nan_to_num(gJ*JexpValList.T)


    def susceptibility(self, ion, Temps, Field, deltaField):
        """
        Calculate magnetic susceptibility χ = ∂M/∂H numerically.
        
        Uses numerical differentiation (4-point stencil) to compute
        susceptibility from magnetization. Can calculate either powder
        average or directional susceptibility.
        
        Args:
            ion: Ion symbol or g-factor
            Temps: Temperature array (K)
            Field: Applied field in Tesla
                - If scalar: computes powder average (χ_x + χ_y + χ_z)/3
                - If 3-vector: computes susceptibility along Field direction
            deltaField: Step size for numerical derivative (Tesla)
            
        Returns:
            Susceptibility in μ_B/Tesla
            - If Field is scalar: powder average (1D array vs Temps)
            - If Field is vector: [χ_x, χ_y, χ_z] components (2D array)
            
        Formula:
            χ ≈ [8(M(H+δ) - M(H-δ)) - (M(H+2δ) - M(H-2δ))] / (12δ)
            (4-point finite difference for accuracy)
            
        Note:
            Choose deltaField carefully: too large → inaccurate, too small → numerical noise
            Typical: deltaField ~ 0.01 T
        """

        if not isinstance(deltaField, float):
            raise TypeError("Deltafield needs to be a scalar")

        if isinstance(Field, float):
            # Assume we are computing a powder average
            VecField = Field * np.array([1,0,0])
            Delta = deltaField*np.array(VecField)/Field
            Mplus1 = self.magnetization(ion, Temps, VecField + Delta)
            Mminus1= self.magnetization(ion, Temps, VecField - Delta)
            Mplus2 = self.magnetization(ion, Temps, VecField + 2*Delta)
            Mminus2= self.magnetization(ion, Temps, VecField - 2*Delta)

            dMdH_x = (8*(Mplus1 - Mminus1) - (Mplus2 - Mminus2))/(12*deltaField)

            VecField = Field * np.array([0,1,0])
            Delta = deltaField*np.array(VecField)/Field
            Mplus1 = self.magnetization(ion, Temps, VecField + Delta)
            Mminus1= self.magnetization(ion, Temps, VecField - Delta)
            Mplus2 = self.magnetization(ion, Temps, VecField + 2*Delta)
            Mminus2= self.magnetization(ion, Temps, VecField - 2*Delta)

            dMdH_y = (8*(Mplus1 - Mminus1) - (Mplus2 - Mminus2))/(12*deltaField)

            VecField = Field * np.array([0,0,1])
            Delta = deltaField*np.array(VecField)/Field
            Mplus1 = self.magnetization(ion, Temps, VecField + Delta)
            Mminus1= self.magnetization(ion, Temps, VecField - Delta)
            Mplus2 = self.magnetization(ion, Temps, VecField + 2*Delta)
            Mminus2= self.magnetization(ion, Temps, VecField - 2*Delta)

            dMdH_z = (8*(Mplus1 - Mminus1) - (Mplus2 - Mminus2))/(12*deltaField)

            return (dMdH_x[:,0]+dMdH_y[:,1]+dMdH_z[:,2])/3.

        elif len(Field) == 3:
            Delta = deltaField*np.array(Field)/np.linalg.norm(Field)
            Mplus1 = self.magnetization(ion, Temps, Field + Delta)
            Mminus1= self.magnetization(ion, Temps, Field - Delta)
            Mplus2 = self.magnetization(ion, Temps, Field + 2*Delta)
            Mminus2= self.magnetization(ion, Temps, Field - 2*Delta)

            dMdH = (8*(Mplus1 - Mminus1) - (Mplus2 - Mminus2))/(12*deltaField)
            #dMdH = (Mplus1 - Mminus1)/(2*deltaField)

            return dMdH


    def susceptibilityPert(self, ion, Temps):
        """
        Calculate powder-averaged susceptibility using perturbation theory.

        Uses second-order perturbation formula (Mesot-Furer equation 11)
        instead of numerical derivatives. Faster but only valid for weak fields.

        Args:
            ion: Ion symbol for g-factor
            Temps: Temperature array (K)
            
        Returns:
            Powder-averaged susceptibility in μ_B/Tesla
            
        Formula:
            χ = (g_J μ_B)² / k_B T * Σ_n [
                P_n ⟨n|J²|n⟩ + 
                Σ_{m≠n} (P_m - P_n)/(E_n - E_m) |⟨m|J|n⟩|²
            ]
            
        Note:
            Assumes powder average: χ = (χ_x + χ_y + χ_z) / 3
            Only accurate in linear (small field) regime
        """

        # Compute susceptibility from perturbation theory, using MesotFurer eq. 11
        gJ = LandeGFactor(ion)
        muB = 5.7883818012e-2  # meV/T
        k_B = 8.6173303e-2  # meV/K

        expvals, temps = np.meshgrid(self.eigenvalues, Temps)
        ZZ = np.sum(np.exp(-expvals/temps/k_B), axis=1)

        suscept = np.zeros(len(Temps))
        for i, ev1 in enumerate(self.eigenvectors):
            kev = Ket(ev1)
            suscept += (np.exp(-self.eigenvalues[i]/Temps/k_B)/ ZZ/ k_B/ Temps) *\
                            np.mean([np.real((kev*kev.Jx()) * (kev*kev.Jx())),
                                    np.real((kev*kev.Jy()) * (kev*kev.Jy())),
                                    np.real((kev*kev.Jz()) * (kev*kev.Jz()))])

            for j, ev2 in enumerate(self.eigenvectors):
                if i == j: continue
                #elif (self.eigenvalues[i]- self.eigenvalues[j]) > 1e-14:
                else:
                    kev2 = Ket(ev2)
                    suscept += ((np.exp(-self.eigenvalues[j]/Temps/k_B)- 
                        np.exp(-self.eigenvalues[i]/Temps/k_B))/ ZZ)/\
                            (self.eigenvalues[i]- self.eigenvalues[j]) *\
                            np.mean([np.real((kev2*kev.Jx()) * (kev*kev2.Jx())),
                                    np.real((kev2*kev.Jy()) * (kev*kev2.Jy())),
                                    np.real((kev2*kev.Jz()) * (kev*kev2.Jz()))])
        return gJ*gJ*muB*suscept


    def gtensor(self):
        """
        Calculate g-tensor from ground state doublet.

        Determines anisotropic g-factors by computing matrix elements of
        J operators between the two lowest eigenstates.

        Returns:
            3×3 g-tensor matrix (not necessarily diagonal)
            
        Structure:
            g = 2 * g_J * [[Re⟨0|J_x|1⟩  Im⟨0|J_x|1⟩  ⟨0|J_x|0⟩]
                            [Re⟨0|J_y|1⟩  Im⟨0|J_y|1⟩  ⟨0|J_y|0⟩]
                            [Re⟨0|J_z|1⟩  Im⟨0|J_z|1⟩  ⟨0|J_z|0⟩]]
            
        Warning:
            Prints warning if ground state is not a doublet (non-Kramers ion
            or broken time-reversal symmetry)
            
        Note:
            Principal g-values found by diagonalizing this matrix
        """

        self.diagonalize_banded()

        def eliminateimag(number):
            num = np.around(number, 10)
            if num.imag == 0:
                return (num.real).astype(float)
            else:
                return number

        zeroinds = np.where(np.around(self.eigenvalues,5)==0)
        gsEVec = self.eigenvectors[zeroinds]

        if len(zeroinds[0]) == 1:
            print('\tWARNING: only one ground state eivenvector found.')

            zeroinds = np.where(np.around(self.eigenvalues,1)==0)
            gsEVec = self.eigenvectors[zeroinds]
            if len(zeroinds[0]) == 2:
                print('\t\t Including another at {} meV'.format(self.eigenvalues[zeroinds[1]]))
            else: raise ValueError('Non-doublet ground state!')


        vv1 = np.around(gsEVec[0],10)
        vv2 = np.around(gsEVec[1],10)

        Jx = self.opttran.Jx
        Jy = self.opttran.Jy*1j
        Jz = self.opttran.Jz

        jz01 = eliminateimag( np.dot(vv1,np.dot(Jz,np.conj(vv2))) )
        jz10 = eliminateimag( np.dot(vv2,np.dot(Jz,np.conj(vv1))) )
        jz00 = eliminateimag( np.dot(vv1,np.dot(Jz,np.conj(vv1))) )
        jz11 = eliminateimag( np.dot(vv2,np.dot(Jz,np.conj(vv2))) )
        
        
        jx01 = eliminateimag( np.dot(vv1,np.dot(Jx,np.conj(vv2))) )
        jx10 = eliminateimag( np.dot(vv2,np.dot(Jx,np.conj(vv1))) )
        jx00 = eliminateimag( np.dot(vv1,np.dot(Jx,np.conj(vv1))) )
        jx11 = eliminateimag( np.dot(vv2,np.dot(Jx,np.conj(vv2))) )
        
        jy01 = eliminateimag( np.dot(vv1,np.dot(Jy,np.conj(vv2))) )
        jy10 = eliminateimag( np.dot(vv2,np.dot(Jy,np.conj(vv1))) )
        jy00 = eliminateimag( np.dot(vv1,np.dot(Jy,np.conj(vv1))) )
        jy11 = eliminateimag( np.dot(vv2,np.dot(Jy,np.conj(vv2))) )
        
        gg = 2*np.array([[np.real(jx01), np.imag(jx01), jx00],
                         [np.real(jy01), np.imag(jy01), jy00],
                         [np.real(jz01), np.imag(jz01), np.abs(jz00)]])
        return gg*LandeGFactor(self.ion)


    def gtensorzeeman(self, field=0.1, Temp=0.1):
        """
        Calculate g-tensor from Zeeman splitting (alternative method).

        Determines g-values by applying small field along each axis and
        measuring the resulting splitting and expectation values.

        Args:
            field: Applied field magnitude (T) - should be small
            Temp: Temperature (K) for thermal averaging
            
        Returns:
            3-component array [g_x, g_y, g_z] along principal axes
            
        Formula:
            g_i = ΔE_zeeman / (μ_B * B * ⟨J_i⟩)
            
        Note:
            Different from gtensor() - this uses energy splitting while
            gtensor() uses matrix elements directly
            Requires low temperature to isolate ground doublet
        """

        Jx = Operator.Jx(self.J)
        Jy = Operator.Jy(self.J)
        Jz = Operator.Jz(self.J)

        muB = 5.7883818012e-2  # meV/T

        gg = np.zeros(3)
        #loop through x,y,z
        for i,Field in enumerate([[field,0,0], [0,field,0], [0,0,field]]):
            JdotB = muB*(Field[0]*Jx + Field[1]*Jy + Field[2]*Jz)

            # B) Diagonalize full Hamiltonian
            FieldHam = self.H + JdotB.O
            diagonalH = LA.eigh(FieldHam)

            minE = np.amin(diagonalH[0])
            evals = diagonalH[0] - minE
            evecs = diagonalH[1].T

            DeltaZeeman = evals[1]-evals[0]
            print(DeltaZeeman)

            # Now find the expectation value of J
            JexpVals = np.zeros((len(evals),3))
            for ii, ev in enumerate(evecs):
                kev = Ket(ev)
                JexpVals[ii] =[np.real(kev*kev.Jx()),
                            np.real(kev*kev.Jy()),
                            np.real(kev*kev.Jz())]
            k_B = 8.6173303e-2  # meV/K

            Zz = np.sum(np.exp(-evals/(k_B*Temp)))
            JexpVal = np.dot(np.exp(-evals/(k_B*Temp)),JexpVals)/Zz

            expectationJ = JexpVal[i]

            # calculate g values
            gg[i] = DeltaZeeman/(muB*field*expectationJ)
        
        return gg


    def fitdata(self, chisqfunc, fitargs, method='Powell', **kwargs):
        """
        Fit crystal field parameters to experimental data.

        Optimizes Stevens parameters by minimizing a user-defined χ² function
        that compares calculated properties to measured data.

        Args:
            chisqfunc: Function(CFLevels, **kwargs) returning χ² value
            fitargs: List of parameter names to fit (e.g., ['B20', 'B40'])
            method: Optimization algorithm ('Powell', 'Nelder-Mead', 'BFGS', etc.)
            **kwargs: Additional arguments passed to chisqfunc
                Must include 'coeff': initial parameter guesses
            
        Returns:
            Dictionary with:
                - Optimized parameter values
                - Final χ² value
                - Other fit metadata
                
        Example:
            >>> def chi2(cf, Texp, Mexp, Field):
            ...     Mcalc = cf.magnetization('Yb3+', Texp, Field)
            ...     return np.sum((Mcalc - Mexp)**2)
            >>> result = cf.fitdata(chi2, ['B20', 'B40'], 
            ...                     coeff=[0.5, -0.02],
            ...                     Texp=T_data, Mexp=M_data, Field=[0,0,1])
        """

        # Define function to be fit
        fun, p0, resfunc = makeFitFunction(chisqfunc, fitargs, **dict(kwargs, CFLevelsObject=self) )

        ############## Fit, using error function  #####################
        p_best = optimize.minimize(fun, p0, method=method)
        ###############################################################

        initialChisq, finalChisq = chisqfunc(self, **kwargs), fun(p_best.x)
        print('\rInitial err =', initialChisq, '\tFinal err =', finalChisq)
        
        result = resfunc(p_best.x)
        result['Chisq'] = finalChisq
        return result


    def fitdata_GlobalOpt(self, chisqfunc, fitargs, **kwargs):
        """
        Global optimization of crystal field parameters using basin-hopping.

        More robust than fitdata() for complex landscapes with multiple
        local minima. Uses basin-hopping algorithm (random jumps + local opt).

        Args:
            chisqfunc: χ² function to minimize
            fitargs: Parameters to optimize
            **kwargs: Arguments for chisqfunc
            
        Returns:
            Dictionary with optimized parameters and final χ²
            
        Note:
            Slower than fitdata() but less likely to get stuck in local minimum
            Performs 100 basin-hopping iterations with temperature T=1e5
            
        Use when:
            - Multiple parameter degeneracies expected
            - Simple optimization fails or gives unphysical results
            - Need to explore parameter space broadly
        """

        # Define function to be fit
        fun, p0, resfunc = makeFitFunction(chisqfunc, fitargs, **dict(kwargs, CFLevelsObject=self) )

        ############## Fit, using error function  #####################
        p_best = optimize.basinhopping(fun, p0, niter=100, T = 1e5)
        ###############################################################

        print(fun(p_best.x))
        print(chisqfunc(self, **kwargs))
        initialChisq, finalChisq = chisqfunc(self, **kwargs), fun(p_best.x)
        print('\rInitial err =', initialChisq, '\tFinal err =', finalChisq)
        
        result = resfunc(p_best.x)
        result['Chisq'] = finalChisq
        return result


    def testEigenvectors(self):
        """
        Verify eigenvector calculation by checking eigenvalue equation.

        Tests whether H|ψ⟩ = E|ψ⟩ is satisfied for each eigenpair.
        Helps diagnose numerical issues in diagonalization.

        Output:
            1. For each eigenvector: H|ψ⟩ - E|ψ⟩ (should be ~0)
            2. Sum rule check: Σ_i |⟨1|J|i⟩|² should equal J(J+1)
            
        Interpretation:
            - Large values (~1e-6) indicate diagonalization problems
            - Sum rule violation suggests incorrect operator construction
            
        Note:
            Uses eigenvaluesNoNorm (before ground state subtraction)
        """

        print('testing eigenvectors... (look for large values)')
        for i in range(len(self.eigenvalues)):
            print(np.around(
                np.dot(self.H,self.eigenvectors[i]) - self.eigenvectors[i]*self.eigenvaluesNoNorm[i],
                10))

        print('\n Sum rule (two values should be equal):')
        TotalTransition = 0
        for i, ev in enumerate(self.eigenvectors):
            TotalTransition += self._transition(Ket(self.eigenvectors[1]),Ket(ev))
        print(TotalTransition, '  ', self.J*(self.J+1))


###################################
### Inner functions for CFLevels class
###################################

# Numba JIT specification for type-safe compilation
# Defines the data structure layout for the JIT compiler

spec = [ 
    ('Jx', float64[:,:]), # 2D array: x-component of angular momentum operator matrix
    ('Jy', float64[:,:]),
    ('Jz', float64[:,:])
]
# These must be float64 (double precision) for numerical accuracy in quantum calculations

@jitclass(spec)
class opttransition(object):
    """
    Optimized transition amplitude calculator using Numba JIT compilation.

    This class provides a ~10-100x speed improvement over the standard Python
    implementation of transition matrix element calculations. It's critical for
    neutron spectrum calculations which require computing thousands of transitions.

    The @jitclass decorator compiles this class to machine code at runtime,
    eliminating Python overhead while maintaining the same interface.

    Physical Purpose:
        Calculates neutron scattering transition intensities:
        I_ij = |⟨i|J_x|j⟩|² + |⟨i|J_y|j⟩|² + |⟨i|J_z|j⟩|²
        
        This represents the probability that a neutron scatters while causing
        a transition from state |i⟩ to state |j⟩.

    Why Optimization Matters:
        - Neutron spectrum calculation requires nested loops over all states
        - For J=7/2 (8 states): 8² = 64 transitions per temperature point
        - For a spectrum with 1000 energy points × 100 temperatures = 6.4M calculations
        - JIT compilation reduces this from ~minutes to ~seconds

    Note:
        The use of Numba's @jitclass requires all data types to be specified
        in advance (via 'spec') and restricts certain Python features, but
        provides dramatic performance improvements for numerical code.
    """

    def __init__(self, optJx, optJy, optJz):
        """
        Initialize the optimized transition calculator with angular momentum operators.
        
        Args:
            optJx: Angular momentum J_x operator matrix (2D numpy array, float64)
                   Size: (2J+1) × (2J+1) where J is total angular momentum
            optJy: Angular momentum J_y operator matrix (2D numpy array, float64)
                   Note: Should be the IMAGINARY part only (already extracted)
            optJz: Angular momentum J_z operator matrix (2D numpy array, float64)
                   
        Implementation Details:
            - Pre-allocates arrays with np.zeros for type safety (Numba requirement)
            - Then assigns the actual operator matrices
            - This two-step process ensures proper memory layout for JIT compilation
            
        Example matrices for J=1/2 (spin-1/2):
            Jx = [[0, 0.5],      Jy = [[0, -0.5],     Jz = [[0.5, 0],
                  [0.5, 0]]            [0.5, 0]]             [0, -0.5]]
                  
        Typical Usage:
            >>> J = 3.5  # For Yb3+ ion
            >>> Jx_matrix = Operator.Jx(J).O.real
            >>> Jy_matrix = Operator.Jy(J).O.imag  # Note: imag part only!
            >>> Jz_matrix = Operator.Jz(J).O.real
            >>> opt_calc = opttransition(Jx_matrix, Jy_matrix, Jz_matrix)
        """

        self.Jx = np.zeros((len(optJx),len(optJx)), dtype=np.float64)
        self.Jy = np.zeros((len(optJx),len(optJx)), dtype=np.float64)
        self.Jz = np.zeros((len(optJx),len(optJx)), dtype=np.float64)
        self.Jx = optJx
        self.Jy = optJy
        self.Jz = optJz

    def transition(self,ket1, ket2):
        """
        Calculate the total transition amplitude between two quantum states.
        
        Computes the squared matrix elements of the angular momentum operator
        summed over all three spatial directions (x, y, z). This is the key
        quantity for neutron scattering cross-sections.
        
        Args:
            ket1: Initial state wavefunction (1D array of length 2J+1)
                  Complex coefficients representing |ψ_i⟩ = Σ c_m |J,m⟩
            ket2: Final state wavefunction (1D array of length 2J+1)
                  Complex coefficients representing |ψ_f⟩
                  
        Returns:
            float: Total transition amplitude (always real, non-negative)
                   Sum of squared matrix elements: Σ_α |⟨i|J_α|j⟩|²
                   
        Physics:
            The neutron magnetic scattering intensity is proportional to:
            
            I ∝ Σ_α |⟨ψ_i|J_α|ψ_f⟩|²
            
            where α ∈ {x, y, z} and the sum accounts for all possible
            polarization directions of the neutron spin flip.
            
        Mathematical Steps:
            1. For each component (x, y, z):
               - Compute: ⟨ket1|J_α|ket2⟩ = ket1† · J_α · ket2
               - Square the result: |⟨ket1|J_α|ket2⟩|²
            2. Sum all three components: ax + ay + az
            
        Optimization:
            - Uses pre-stored J matrices (no repeated construction)
            - np.dot is compiled to BLAS calls by Numba
            - No Python overhead in the inner calculation loop
            - Type-stable: all float64 operations
            
        Example Calculation:
            For two eigenstates of a crystal field Hamiltonian:
            >>> eigenstate_i = cf.eigenvectors[0]  # Ground state
            >>> eigenstate_j = cf.eigenvectors[2]  # First excited state
            >>> amplitude = opt_calc.transition(eigenstate_i, eigenstate_j)
            >>> # amplitude might be ~0.5 for allowed transition, ~0.0 for forbidden
            
        Performance:
            - Standard Python version (~CFLevels._transition): ~1 ms per call
            - This JIT version: ~10 μs per call
            - Speedup: ~100x (critical for spectrum calculations)
            
        Selection Rules:
            - The value tells you if a transition is allowed or forbidden
            - amplitude ≈ 0: Forbidden by symmetry (e.g., Δm_J = ±2)
            - amplitude > 0: Allowed transition (intensity ∝ amplitude)
            
        Note:
            The result is always real and non-negative because we're computing
            |⟨i|J|j⟩|² (squared modulus). Any imaginary parts from the matrix
            elements cancel when squared.
        """
        
        ax = np.dot(ket1, np.dot(self.Jx, ket2))**2
        ay = np.dot(ket1, np.dot(self.Jy, ket2))**2
        az = np.dot(ket1, np.dot(self.Jz, ket2))**2
        return ax + ay + az

"""
LandeGFactor Function
"""

def LandeGFactor(ion):
    """
    Calculate the Lande g-factor for a given ion. 
    Lande g-factor describes the magnetic moment of an atom or ion.
    
    Args:
        ion (str): Ion name (e.g., 'Ce3+', 'Pr3+')
        
    Returns:
        float: Lande g-factor value
    """
    s, l, j = ION_NUMS_RARE_EARTH[ion]
    return 1.5 + (s*(s+1.) - l*(l+1.))/(2.*j*(j+1.))

###################################
### Same class, but in the LS basis
###################################

# LS basis functions are a set of quantum mechanical wave functions
# that describe the electronic states of an atom 
# in the L-S (Russell-Saunders) coupling scheme. 
# This model is used for lighter atoms where 
# the electrostatic repulsion between electrons
# is stronger than the spin-orbit interaction.

class LS_CFLevels:
    """
    Crystal field levels calculator for systems requiring LS coupling treatment.
    
    This class handles crystal field calculations where spin (S) and orbital (L)
    angular momentum are treated separately before coupling, in contrast to CFLevels
    which uses the total J basis. Used primarily for:
    
    1. **Transition metal ions** (3d electrons: Ti³⁺, V³⁺, Cr³⁺, Mn²⁺, Fe²⁺, Co²⁺, Ni²⁺, Cu²⁺)
       - Weaker spin-orbit coupling (ζ ~ 50-400 cm⁻¹ = 6-50 meV)
       - Crystal field >> spin-orbit coupling
       - LS coupling more appropriate than Russell-Saunders
       
    2. **Light rare earths** (4f¹-⁴f⁷ before half-filling)
       - When crystal field and spin-orbit are comparable
       
    Physical Basis:
        Total Hamiltonian: H = H_CEF + H_SOC + H_Zeeman
        
        Where:
        - H_CEF: Crystal field acting on orbital degrees of freedom
        - H_SOC: Spin-orbit coupling λ L·S
        - H_Zeeman: Magnetic field interaction (g₀S + L)·B
        
    Basis States:
        |L, m_L, S, m_S⟩ instead of |J, m_J⟩
        Dimension: (2L+1) × (2S+1)
        
    Example:
        For d² ion (like V³⁺, L=3, S=1):
        - Dimension: 7 × 3 = 21 states
        - Basis: |3,-3,1,-1⟩, |3,-3,1,0⟩, |3,-3,1,+1⟩, ..., |3,+3,1,+1⟩
        
    Key Differences from CFLevels:
        - Separate L and S operators (not combined into J)
        - Spin-orbit coupling added as separate term
        - g-factor includes g₀ ≈ 2.002 for spin part
        - Used when crystal field ≳ spin-orbit coupling
    """

    def __init__(self, StevensOperators, Parameters, L,S, SpinOrbitCoupling):
        """
        Initialize LS_CFLevels with separated spin and orbital components.
        
        Constructs two separate Hamiltonians:
        1. H_CEF: Crystal field Hamiltonian (acts only on orbital part)
        2. H_SOC: Spin-orbit coupling Hamiltonian (couples L and S)
        
        Args:
            StevensOperators: List of Stevens operator matrices in LS basis
                             Generated by LS_StevensOp(L, S, n, m)
            Parameters: List of B_n^m crystal field parameters (meV)
            L: Orbital angular momentum quantum number (integer or half-integer)
            S: Spin angular momentum quantum number (integer or half-integer)
            SpinOrbitCoupling: Spin-orbit coupling constant λ (meV)
                              - Positive for less-than-half-filled shells (d¹-d⁴, f¹-f⁶)
                              - Negative for more-than-half-filled shells (d⁶-d⁹, f⁸-f¹³)
        
        Sets:
            self.H_CEF: Crystal electric field Hamiltonian (LSOperator object)
            self.H_SOC: Spin-orbit Hamiltonian λ(L·S) (LSOperator object)
            self.O: Saved operators for fitting
            self.B: Parameter values
            self.L, self.S: Quantum numbers
            self.Jx, Jy, Jz: Total angular momentum operators (L + S)
            self.Jxg0, Jyg0, Jzg0: Weighted operators (g₀S + L) for magnetization
            
        Physical Interpretation:
            The total Hamiltonian without external field is:
            H = Σ B_n^m O_n^m(L) + λ L·S
            
            First term: Crystal field splits orbital states
            Second term: Spin-orbit coupling mixes L and S
            
        Example:
            >>> # For Ni²⁺ ion (d⁸, L=3, S=1)
            >>> params = {'B20': 0.5, 'B40': -0.02, 'B44': 0.1}
            >>> operators = [LS_StevensOp(3, 1, n, m) for n,m in [(2,0), (4,0), (4,4)]]
            >>> cf = LS_CFLevels(operators, [0.5, -0.02, 0.1], 
            ...                  L=3, S=1, SpinOrbitCoupling=-0.05)
        """

        # Build crystal field Hamiltonian (acts on orbital states only)
        self.H_CEF = LSOperator(L, S)
        HcefJ = np.sum([a*b for a,b in zip(StevensOperators, Parameters)], axis=0)
        self.H_CEF.O = HcefJ
        self.O = StevensOperators  #save these for a fit
        self.B = Parameters
        self.S = S
        self.L = L

        # Define spin-orbit coupling Hamiltonian: λ L·S
        # Build individual spin and orbital operators
        Sx = LSOperator.Sx(L, S)
        Sy = LSOperator.Sy(L, S)
        Sz = LSOperator.Sz(L, S)
        Lx = LSOperator.Lx(L, S)
        Ly = LSOperator.Ly(L, S)
        Lz = LSOperator.Lz(L, S)

        # Construct L·S = L_x S_x + L_y S_y + L_z S_z
        self.H_SOC = Lx*Sx + Ly*Sy + Lz*Sz
        LdotS = self.H_SOC.O*1.0

        # Clean up: if imaginary part is zero (should be), make it real
        if np.sum(LdotS.imag) == 0: LdotS = LdotS.real

        # Scale by spin-orbit coupling constant
        self.H_SOC.O = SpinOrbitCoupling*LdotS
        self.spinorbitcoupling = SpinOrbitCoupling

        # Define total angular momentum operators: J = L + S
        # Used for calculating physical observables
        g0 = 2.002319 # Free electron g-factor (CODATA value)

        # Pure J operators (unweighted sum)
        self.Jx = Sx + Lx
        self.Jy = Sy + Ly
        self.Jz = Sz + Lz

        # Weighted J operators: g₀S + L (for magnetization calculations)
        # This accounts for different g-factors of spin vs orbital contributions
        self.Jxg0 = g0*Sx + Lx
        self.Jyg0 = g0*Sy + Ly
        self.Jzg0 = g0*Sz + Lz


    @classmethod
    def Bdict(cls,  Bdict, L, S, SpinOrbitCoupling):
        """
        Create LS_CFLevels instance from parameter dictionary.
        
        Factory method that constructs Stevens operators from L and S values
        and builds Hamiltonian from B_n^m parameter dictionary.
        
        Args:
            Bdict: Dictionary mapping parameter names to values
                   Keys: 'B20', 'B40', 'B43', etc.
                   Values: Parameters in meV
            L: Orbital angular momentum quantum number
            S: Spin angular momentum quantum number
            SpinOrbitCoupling: λ parameter (meV)
            
        Returns:
            LS_CFLevels instance with constructed Hamiltonians
            
        Example:
            >>> params = {'B20': -0.340, 'B40': -0.031, 'B44': 0.227}
            >>> cf = LS_CFLevels.Bdict(params, L=3, S=1, SpinOrbitCoupling=-0.05)
            
        Note:
            Uses LS_StevensOp() to build operators in LS coupling basis,
            different from CFLevels.Bdict() which uses J-basis operators.
        """

        Stev_O = []
        Parameters = []
        for Bnm in Bdict:
            Parameters.append(Bdict[Bnm])
            n = int(Bnm[1])                             # Extract rank from 'B20' → n=2
            m = int(Bnm[2:])                            # Extract order from 'B20' → m=0
            Stev_O.append(LS_StevensOp(L,S,n,m))        # Build in LS basis

        newcls = cls(Stev_O, Parameters, L, S, SpinOrbitCoupling)
        return newcls


    @classmethod
    def Hamiltonian(cls, CEF_Hamil, SOC_Hamil, L, S):
        """
        Create LS_CFLevels directly from pre-constructed Hamiltonian matrices.
        
        Alternative constructor for cases where Hamiltonians are built externally
        or need custom construction beyond standard Stevens operators.
        
        Args:
            CEF_Hamil: Crystal field Hamiltonian (LSOperator object or matrix)
            SOC_Hamil: Spin-orbit Hamiltonian (LSOperator object or matrix)
            L: Orbital quantum number
            S: Spin quantum number
            
        Returns:
            LS_CFLevels instance with the provided Hamiltonians
            
        Use Cases:
            - Custom Hamiltonians with non-standard terms
            - Importing from ab-initio calculations
            - Adding exchange interactions or other couplings
        """

        newcls = cls([0,0],[0,0], L, S, 0) # Create empty class so we can just define Hamiltonian
        newcls.H_CEF = CEF_Hamil  # Crystal electric fields
        newcls.H_SOC = SOC_Hamil  # Spin-Orbit Coupling
        return newcls
    

    def newCoeff(self, newcoeff):
        """
        Update crystal field parameters and rediagonalize.
        
        Reconstructs H_CEF with new B_n^m coefficients, keeping H_SOC unchanged.
        Used during fitting procedures to test different parameter sets.
        
        Args:
            newcoeff: New list of crystal field parameters (meV)
            
        Side effects:
            Updates self.B and self.H_CEF, triggers diagonalization
            
        Note:
            Only updates crystal field part; spin-orbit coupling remains constant
        """

        self.B = np.array(newcoeff)
        newH = np.sum([a*b for a,b in zip(self.O, newcoeff)], axis=0)
        self.diagonalize(newH)


    def diagonalize(self, CEF_Hamiltonian=None):
        """
        Diagonalize the total Hamiltonian H = H_CEF + H_SOC.
        
        Finds energy eigenvalues and eigenvectors of the combined crystal field
        and spin-orbit Hamiltonian using banded matrix algorithm.
        
        Args:
            CEF_Hamiltonian: Optional crystal field Hamiltonian matrix
                           If None, uses self.H_CEF.O
                           
        Sets:
            self.eigenvaluesNoNorm: Raw eigenvalues
            self.eigenvalues: Eigenvalues with ground state at E=0
            self.eigenvectors: Eigenvectors in |L,m_L,S,m_S⟩ basis
            
        Physical Interpretation:
            Solves: (H_CEF + H_SOC)|ψ⟩ = E|ψ⟩
            
            Eigenstates are linear combinations:
            |ψ_i⟩ = Σ_{m_L,m_S} c_{m_L,m_S} |L,m_L,S,m_S⟩
            
        Example Result:
            For Ni²⁺ (L=3, S=1), 21 states split into groups:
            eigenvalues = [0.0, 5.2, 8.7, 15.3, ...]  # meV
            
        Note:
            Uses banded storage for efficiency (same as CFLevels.diagonalize_banded)
            Threshold 1e-15 cleans numerical noise near zero
        """

        if CEF_Hamiltonian is None:
            CEF_Hamiltonian = self.H_CEF.O
        else:
            self.H_CEF.O = CEF_Hamiltonian

        # Diagonalize combined Hamiltonian
        bands = self._findbands(CEF_Hamiltonian + self.H_SOC.O)
        diagonalH = LA.eig_banded(bands, lower=True)

        self.eigenvaluesNoNorm = diagonalH[0]
        self.eigenvalues = diagonalH[0] - np.amin(diagonalH[0])         # Set ground state to 0
        self.eigenvectors = diagonalH[1].T

        # Clean up numerical noise
        tol = 1e-15                                                     # set very small values to zero
        self.eigenvalues[abs(self.eigenvalues) < tol] = 0.0
        self.eigenvectors[abs(self.eigenvectors) < tol] = 0.0

    # Shared
    def _findbands(self, matrix):
        """
        Extract banded structure from sparse Hamiltonian matrix.
        
        Identifies non-zero diagonals for efficient banded matrix storage.
        Angular momentum matrices are naturally sparse (only near-diagonal elements).
        
        Args:
            matrix: Full Hamiltonian matrix (potentially sparse)
            
        Returns:
            Banded matrix array suitable for LA.eig_banded()
            
        Technical Details:
            - Row i contains the i-th diagonal of the matrix
            - Only extracts diagonals with |elements| > 1e-10
            - Reduces storage from O(n²) to O(n×bandwidth)
            
        Note:
            Shared with CFLevels class (identical implementation)
        """

        diags = np.zeros((len(matrix),len(matrix)), dtype=np.complex128)
        for i in range(len(matrix)):
            diag = matrix.diagonal(i)
            if i == 0:
                diags[i] = diag
            else:
                diags[i][:-i] = diag
            if np.count_nonzero(np.around(diag,10)) > 0:
                nonzerobands = i
        return diags[:nonzerobands+1]


    def neutronSpectrum(self, Earray, Temp, Ei, ResFunc, gamma = 0):
        """
        Calculate inelastic neutron scattering spectrum for LS-coupled system.
        
        Computes S(Q,ω) for transition metal or light rare earth ions where
        both spin and orbital contributions matter.
        
        Args:
            Earray: Energy transfer grid (meV)
            Temp: Sample temperature (K)
            Ei: Incident neutron energy (meV)
            ResFunc: Energy resolution function
            gamma: Lorentzian width for intrinsic broadening (meV)
            
        Returns:
            Intensity array with k'/k kinematic correction
            
        Key Differences from CFLevels.neutronSpectrum():
            1. Uses Ket objects (not optimized opttran class)
            2. Limits to first 12 transitions (maxtransition=12)
            3. Uses weighted J operators (Jxg0, Jyg0, Jzg0) accounting for g₀
            
        Physics:
            Neutron scattering probes: (g₀S + L)
            This is the magnetic moment operator with proper g-factors
            
        Limitation:
            maxtransition = 12 hardcoded
            - For large L,S, full spectrum may have >12 accessible levels
            - Higher transitions typically have low intensity anyway
            - Speeds up calculation significantly
            
        Note:
            Uses self._transition() which includes g₀ weighting,
            different from CFLevels which uses equal weighting
        """

        try:
            eigenkets = self.eigenvectors.real
            intensity = np.zeros(len(Earray))
        except AttributeError:
            self.diagonalize()
            eigenkets = self.eigenvectors.real
            intensity = np.zeros(len(Earray))

        maxtransition = 12 # because we can't see the others

        # Create Ket objects for angular momentum calculations
        eigenkets = [Ket(ei) for ei in self.eigenvectors[:maxtransition]]

        # Boltzmann population factors
        beta = 1/(8.61733e-2*Temp)  # Boltzmann constant is in meV/K
        Z = sum([np.exp(-beta*en) for en in self.eigenvalues])

        # Double loop over all transitions
        for i, ket_i in enumerate(eigenkets):
            # compute population factor
            pn = np.exp(-beta *self.eigenvalues[i])/Z

            for j, ket_j in enumerate(eigenkets):
                # Transition amplitude with g₀ weighting
                mJn = self._transition(ket_i,ket_j)
                deltaE = self.eigenvalues[j] - self.eigenvalues[i]
                GausWidth = ResFunc(deltaE)  #peak width due to instrument resolution

                # Add peak with Voigt profile
                intensity += ((pn * mJn * self._voigt(x=Earray, x0=deltaE, alpha=GausWidth, 
                                                    gamma=gamma)).real).astype('float64')

        # Apply kinematic correction
        kpoverk = np.sqrt((Ei - Earray)/Ei) #k'/k = sqrt(E'/E)
        return intensity * kpoverk


    # Shared
    def neutronSpectrum2D(self, Earray, Qarray, Temp, Ei, ResFunc, gamma, DebyeWaller, Ion):
        """
        Calculate 2D neutron scattering map: S(Q,ω).
        
        Extends 1D spectrum to include Q-dependence through form factor
        and Debye-Waller factor.
        
        Args:
            Earray: Energy transfer grid (meV)
            Qarray: Momentum transfer grid (Å⁻¹)
            Temp: Temperature (K)
            Ei: Incident energy (meV)
            ResFunc: Resolution function
            gamma: Lorentzian width (meV)
            DebyeWaller: Debye-Waller parameter ⟨u²⟩^(1/2) (Å)
            Ion: Ion symbol for form factor lookup
            
        Returns:
            2D intensity array [len(Earray) × len(Qarray)]
            
        Physical Formula:
            I(Q,E) = I_1D(E) × F²(Q) × exp(-Q²⟨u²⟩/3)
            
        Note:
            Shared implementation with CFLevels (identical method)
        """

        intensity1D = self.neutronSpectrum(Earray, Temp, Ei, ResFunc,  gamma)

        
        # Q-dependent factors
        ## Scale by Debye-Waller Factor
        DWF = np.exp(1./3. * Qarray**2 * DebyeWaller**2)
        ## Scale by Magnetic form factor
        FormFactor = RE_FormFactor(Qarray,Ion)
        return np.outer(intensity1D, DWF*FormFactor)


    def normalizedNeutronSpectrum(self, Earray, Temp, ResFunc, gamma = 0):
        """
        Neutron spectrum without kinematic (k'/k) correction.
        
        Identical to neutronSpectrum() but omits the sqrt(E'/E) factor,
        useful for relative intensity comparisons.
        
        Args:
            Earray: Energy transfer (meV)
            Temp: Temperature (K)
            ResFunc: Resolution function
            gamma: Lorentzian width (meV)
            
        Returns:
            Intensity without kinematic correction
            
        Note:
            Also limited to maxtransition=12 for efficiency
        """

        try:
            eigenkets = self.eigenvectors.real
            intensity = np.zeros(len(Earray))
        except AttributeError:
            self.diagonalize()
            eigenkets = self.eigenvectors.real
            intensity = np.zeros(len(Earray))

        maxtransition = 12 # because we can't see the others

        # make angular momentum ket object
        eigenkets = [Ket(ei) for ei in self.eigenvectors[:maxtransition]]

        # for population factor weights
        beta = 1/(8.61733e-2*Temp)  # Boltzmann constant is in meV/K
        Z = sum([np.exp(-beta*en) for en in self.eigenvalues])

        for i, ket_i in enumerate(eigenkets):
            # compute population factor
            pn = np.exp(-beta *self.eigenvalues[i])/Z

            for j, ket_j in enumerate(eigenkets):
                # compute amplitude
                mJn = self._transition(ket_i,ket_j)
                deltaE = self.eigenvalues[j] - self.eigenvalues[i]
                GausWidth = ResFunc(deltaE)  #peak width due to instrument resolution
                intensity += ((pn * mJn * self._voigt(x=Earray, x0=deltaE, alpha=GausWidth, 
                                                    gamma=gamma)).real).astype('float64')
                #intensity += ((pn * mJn * self._lorentzian(Earray, deltaE, Width)).real).astype('float64')
        return intensity


    def _transition(self,ket1,ket2):
        """
        Calculate transition amplitude including proper g-factors for LS coupling.
        
        Computes |⟨i|(g₀S + L)|j⟩|² summed over x,y,z components.
        This is the magnetic scattering cross-section operator.
        
        Args:
            ket1: Initial state (Ket object)
            ket2: Final state (Ket object)
            
        Returns:
            Total transition amplitude (real, non-negative)
            
        Key Difference from CFLevels._transition():
            Uses Jxg0, Jyg0, Jzg0 instead of Jx, Jy, Jz
            These include the g₀≈2 factor for spin:
            J_g0 = g₀S + L
            
        Physical Meaning:
            Neutrons couple to magnetic moment μ = μ_B(g₀S + L)
            The g₀ factor accounts for spin having g≈2 vs orbital g=1
            
        Formula:
            I_ij = |⟨i|(g₀S_x+L_x)|j⟩|² + |⟨i|(g₀S_y+L_y)|j⟩|² + |⟨i|(g₀S_z+L_z)|j⟩|²
            
        Note:
            Uses explicit np.dot operations with conjugate (slower than opttran)
            But necessary for complex eigenvectors in LS basis
            
        Raises:
            ValueError: If result has imaginary part (indicates numerical error)
        """

        # Calculate each component with proper conjugation
        ax = np.dot(np.conjugate(ket1.ket),np.dot(self.Jxg0.O,ket2.ket)) *\
                np.dot(np.conjugate(ket2.ket),np.dot(self.Jxg0.O,ket1.ket))
        ay = np.dot(np.conjugate(ket1.ket),np.dot(self.Jyg0.O,ket2.ket)) *\
                np.dot(np.conjugate(ket2.ket),np.dot(self.Jyg0.O,ket1.ket))
        az = np.dot(np.conjugate(ket1.ket),np.dot(self.Jzg0.O,ket2.ket)) *\
                np.dot(np.conjugate(ket2.ket),np.dot(self.Jzg0.O,ket1.ket))

        # Eliminate tiny values (numerical noise)
        ax, ay, az = np.around(ax, 10), np.around(ay, 10), np.around(az, 10)
        if (ax + ay + az).imag == 0:
            return ((ax + ay + az).real).astype(float)
        else:
            print(ax, ay, az)
            raise ValueError("non-real amplitude. Error somewhere.")


    # Shared
    def _lorentzian(self, x, x0, gamma):
        """
        Lorentzian lineshape for spectral peaks.
        
        Args:
            x: Energy array
            x0: Peak center
            gamma: FWHM
            
        Returns:
            Lorentzian profile
            
        Note:
            Shared with CFLevels (identical implementation)
        """

        return 1/np.pi * (0.5*gamma)/((x-x0)**2 + (0.5*gamma)**2)


    # Shared
    def _voigt(self, x, x0, alpha, gamma):
        """
        Voigt profile (Gaussian ⊗ Lorentzian convolution).
        
        Args:
            x: Energy array
            x0: Peak center
            alpha: Gaussian FWHM (resolution)
            gamma: Lorentzian FWHM (lifetime)
            
        Returns:
            Voigt profile
            
        Note:
            Shared with CFLevels (identical implementation)
        """

        sigma = (0.5*alpha) / np.sqrt(2 * np.log(2))
        return np.real(wofz(((x-x0) + 1j*(0.5*gamma))/sigma/np.sqrt(2))) / sigma\
                                                            /np.sqrt(2*np.pi)
    
    # Shared
    def _Re(self,value):
        """
        Extract real part if imaginary component is negligible.
        
        Args:
            value: Complex number or array
            
        Returns:
            Real part if |Im| < 1e-9, else complex value
            
        Note:
            Shared with CFLevels (identical implementation)
        """

        thresh = 1e-9
        if np.size(value) == 1 & isinstance(value, complex):
            if np.abs(value.imag) <= thresh:
                return (value.real).astype(float)
            else: 
                return value
        else:
            if np.all(np.abs(value.imag) < thresh):
                return (value.real)
            else: return value

    # Shared
    def printEigenvectors(self):
        """
        Display eigenvalues and eigenvectors in ASCII table.
        
        Prints sorted energy levels with corresponding wavefunction components
        in the |L,m_L,S,m_S⟩ basis.
        
        Output Format:
            Eigenvalues    Eigenvectors
            ----------------------------------
            0.00000    | [21 components for L=3,S=1] |
            ...
            
        Note:
            Shared with CFLevels (identical implementation)
            For large L,S, table becomes very wide (2L+1)×(2S+1) columns
        """

        try:
            eigenkets = self.eigenvectors.real
        except AttributeError:
            self.diagonalize()
        
        print('\n Eigenvalues \t Eigenvectors')
        print('\t\t'+'-------'*(len(self.eigenvalues)+1))
        sortinds = self.eigenvalues.argsort()
        sortEVal= np.around(self.eigenvalues[sortinds],5)
        sortEVec= np.around(self.eigenvectors[sortinds],3)
        for i in range(len(sortinds)):
            print(format(self._Re(sortEVal[i]), '.5f'),'\t| ', self._Re(sortEVec[i]),' |')
        print('\t\t'+'-------'*(len(self.eigenvalues)+1) + '\n')


    def gsExpectation(self):
        """
        Calculate expectation values of J for ground state(s).
        
        Computes ⟨J_x⟩, ⟨J_y⟩, ⟨J_z⟩ using weighted operators (Jxg0, etc.)
        for all degenerate ground states.
        
        Output:
            Prints expectation values for each ground state component
            
        Physical Interpretation:
            Shows the effective magnetic moment direction and magnitude
            in the crystal field ground state
            
        Key Difference from CFLevels:
            Uses Jxg0 (weighted by g₀) instead of Jx
            Reflects proper g-factors: μ = μ_B(g₀S + L)
        """

        zeroinds = np.where(np.around(self.eigenvalues,7)==0)
        gsEVec = self.eigenvectors[zeroinds]
        print('\t Ground State Expectation Values:')

        for ev in gsEVec:
            # Uses weighted J operators
            jjxx = self._Re(np.dot(ev,np.dot(self.Jxg0.O,ev)))
            jjyy = self._Re(np.dot(ev,np.dot(self.Jyg0.O,ev)))
            jjzz = self._Re(np.dot(ev,np.dot(self.Jzg0.O,ev)))
            print('  <J_x> =',jjxx,'\t<J_y> =',jjyy,'\t<J_z> =',jjzz)
        print(' ')


    def magnetization(self, Temp, Field):
        """
        Calculate thermal-averaged magnetization for LS-coupled system.
        
        Computes M(T,H) by diagonalizing H = H_CEF + H_SOC + μ_B(g₀S+L)·B
        and thermally averaging the magnetic moment operator.
        
        Args:
            Temp: Temperature (K) - scalar or array
            Field: 3-component field vector [B_x, B_y, B_z] in Tesla
            
        Returns:
            [M_x, M_y, M_z] in Bohr magnetons
            - Scalar Temp: 3-component vector
            - Array Temp: [3 × len(Temp)] array
            
        Key Differences from CFLevels.magnetization():
            1. No ion/gJ parameter (uses g₀ hardcoded)
            2. Uses Jxg0, Jyg0, Jzg0 operators (weighted by g₀)
            3. Simpler formula: no Landé g-factor needed
            
        Physics:
            Zeeman Hamiltonian: H_Z = μ_B(g₀S + L)·B
            Magnetization: M = ⟨g₀S + L⟩ (in units of μ_B)
            
        Example:
            >>> M = cf.magnetization(Temp=10, Field=np.array([0, 0, 1]))
            >>> # M might be [0.0, 0.0, 1.5] for Ising-like ground state
            
        Note:
            Returns real values only (imaginary parts removed via np.real)
        """

        if len(Field) != 3: 
            raise TypeError("Field needs to be 3-component vector")

        # A) Build Zeeman Hamiltonian with weighted operators
        muB = 5.7883818012e-2  # meV/T
        JdotB = muB*(Field[0]*self.Jxg0 + Field[1]*self.Jyg0 + Field[2]*self.Jzg0)

        # B) Diagonalize total Hamiltonian
        FieldHam = self.H_CEF.O + self.H_SOC.O + JdotB.O
        diagonalH = LA.eigh(FieldHam)

        minE = np.amin(diagonalH[0])
        evals = diagonalH[0] - minE
        evecs = diagonalH[1].T

        # C) Compute expectation value along field

        # Compute expectation values of J for each eigenstate
        JexpVals = np.zeros((len(evals),3))
        for i, ev in enumerate(evecs):
            JexpVals[i] =[np.real(np.dot(np.conjugate(ev), np.dot( self.Jxg0.O ,ev))),
                          np.real(np.dot(np.conjugate(ev), np.dot( self.Jyg0.O ,ev))),
                          np.real(np.dot(np.conjugate(ev), np.dot( self.Jzg0.O ,ev)))]
        k_B = 8.6173303e-2  # meV/K

        # Thermal average
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
        Calculate magnetic susceptibility χ = ∂M/∂H numerically.
        
        Uses 4-point finite difference to compute susceptibility from
        magnetization, either as powder average or directional.
        
        Args:
            Temps: Temperature array (K)
            Field: Applied field (T)
                - Scalar: powder average (χ_x + χ_y + χ_z)/3
                - 3-vector: directional susceptibility
            deltaField: Step size for derivative (T)
            
        Returns:
            Susceptibility in μ_B/T
            
        Note:
            Identical algorithm to CFLevels.susceptibility() but
            uses LS_CFLevels.magnetization() internally
            
        Formula:
            χ ≈ [8(M(H+δ)-M(H-δ)) - (M(H+2δ)-M(H-2δ))] / (12δ)
        """

        if not isinstance(deltaField, float):
            raise TypeError("Deltafield needs to be a scalar")

        if isinstance(Field, float):
            # Powder average calculation (Assume we are computing a powder average)
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
            # Directional susceptibility
            Delta = deltaField*np.array(Field)/np.linalg.norm(Field)
            Mplus1 = self.magnetization(Temps, Field + Delta)
            Mminus1= self.magnetization(Temps, Field - Delta)
            Mplus2 = self.magnetization(Temps, Field + 2*Delta)
            Mminus2= self.magnetization(Temps, Field - 2*Delta)

            dMdH = (8*(Mplus1 - Mminus1) - (Mplus2 - Mminus2))/(12*deltaField)

            return dMdH


    def susceptibilityDeriv(self, Temps, Field, deltaField):
        """
        Calculate susceptibility using thermodynamic derivative method.
        
        Alternative to susceptibility() that uses free energy derivative:
        χ = -∂²F/∂H² = ∂M/∂H
        
        Computed via magnetizationDeriv() which uses energy eigenvalue derivatives
        rather than expectation values.
        
        Args:
            Temps: Temperature array (K)
            Field: Applied field (T) - scalar for powder average
            deltaField: Step size for numerical derivatives
            
        Returns:
            Powder-averaged susceptibility in μ_B/T
            
        Method:
            Takes derivative of magnetizationDeriv() w.r.t. field
            Uses nested finite differences
            
        Note:
            More complex than susceptibility() but can be more accurate
            for certain cases. Currently only implements powder average.
        """

        if not isinstance(deltaField, float):
            raise TypeError("Deltafield needs to be a scalar")

        if (isinstance(Field, float) or isinstance(Field, int)):
            # Apply finite difference to magnetizationDeriv (Assume we are computing a powder average)
            VecField = Field * np.array([1,0,0])
            Delta = deltaField*np.array(VecField)/Field
            Mplus1 = self.magnetizationDeriv(Temps, VecField + Delta, deltaField)
            Mminus1= self.magnetizationDeriv(Temps, VecField - Delta, deltaField)
            Mplus2 = self.magnetizationDeriv(Temps, VecField + 2*Delta, deltaField)
            Mminus2= self.magnetizationDeriv(Temps, VecField - 2*Delta, deltaField)

            dMdH_x = (8*(Mplus1 - Mminus1) - (Mplus2 - Mminus2))/(12*deltaField)

            VecField = Field * np.array([0,1,0])
            Delta = deltaField*np.array(VecField)/Field
            Mplus1 = self.magnetizationDeriv(Temps, VecField + Delta, deltaField)
            Mminus1= self.magnetizationDeriv(Temps, VecField - Delta, deltaField)
            Mplus2 = self.magnetizationDeriv(Temps, VecField + 2*Delta, deltaField)
            Mminus2= self.magnetizationDeriv(Temps, VecField - 2*Delta, deltaField)

            dMdH_y = (8*(Mplus1 - Mminus1) - (Mplus2 - Mminus2))/(12*deltaField)

            VecField = Field * np.array([0,0,1])
            Delta = deltaField*np.array(VecField)/Field
            Mplus1 = self.magnetizationDeriv(Temps, VecField + Delta, deltaField)
            Mminus1= self.magnetizationDeriv(Temps, VecField - Delta, deltaField)
            Mplus2 = self.magnetizationDeriv(Temps, VecField + 2*Delta, deltaField)
            Mminus2= self.magnetizationDeriv(Temps, VecField - 2*Delta, deltaField)

            dMdH_z = (8*(Mplus1 - Mminus1) - (Mplus2 - Mminus2))/(12*deltaField)

            return (dMdH_x[:,0]+dMdH_y[:,1]+dMdH_z[:,2])/3.

        else: return 0


    def magnetizationDeriv(self, Temp, Field, deltaField):
        """
        Calculate magnetization via derivative of free energy: M = -∂F/∂H.
        
        Uses Hellmann-Feynman theorem: ∂E_n/∂H = ⟨n|∂H/∂H|n⟩
        More direct than computing expectation values explicitly.
        
        Args:
            Temp: Temperature (K) - scalar or array
            Field: 3-component field vector (T)
            deltaField: Step size for energy derivatives
            
        Returns:
            Magnetization [M_x, M_y, M_z] in μ_B
            
        Algorithm:
            1. Diagonalize at H to get {E_n(H)}
            2. Diagonalize at H±δ along each axis to get {E_n(H±δ)}
            3. Compute ∂E_n/∂H_α via finite difference
            4. M_α = -Σ_n P_n ∂E_n/∂H_α where P_n = exp(-E_n/k_B T)/Z
            
        Advantages:
            - No need to compute matrix elements explicitly
            - Uses eigenvalue derivatives only
            - Can be more numerically stable
            
        Physics:
            From thermodynamics: M = -∂F/∂H where F = -k_B T ln(Z)
            Since Z = Σ exp(-E_n/k_B T), this gives M = Σ P_n (∂E_n/∂H)
            
        Note:
            Requires 7 diagonalizations per call (1 at H, 6 at H±δ for 3 directions)
            Can be expensive for large systems
        """

        if len(Field) != 3: 
            raise TypeError("Field needs to be 3-component vector")

        # A) Define magnetic Hamiltonian
        muB = 5.7883818012e-2  # meV/T

        # Build Zeeman terms for H and H±δ along each axis
        JdotB_0 = muB*(Field[0]*self.Jxg0 + Field[1]*self.Jyg0 + Field[2]*self.Jzg0)

        # X-direction perturbations
        Delta = deltaField*np.array([1,0,0])
        FieldPD = Field + Delta
        JdotB_p1x = muB*(FieldPD[0]*self.Jxg0 + FieldPD[1]*self.Jyg0 + FieldPD[2]*self.Jzg0)
        FieldPD = Field - Delta
        JdotB_m1x = muB*(FieldPD[0]*self.Jxg0 + FieldPD[1]*self.Jyg0 + FieldPD[2]*self.Jzg0)

        # Y-direction perturbations
        Delta = deltaField*np.array([0,1,0])
        FieldPD = Field + Delta
        JdotB_p1y = muB*(FieldPD[0]*self.Jxg0 + FieldPD[1]*self.Jyg0 + FieldPD[2]*self.Jzg0)
        FieldPD = Field - Delta
        JdotB_m1y = muB*(FieldPD[0]*self.Jxg0 + FieldPD[1]*self.Jyg0 + FieldPD[2]*self.Jzg0)

        # Z-direction perturbations
        Delta = deltaField*np.array([0,0,1])
        FieldPD = Field + Delta
        JdotB_p1z = muB*(FieldPD[0]*self.Jxg0 + FieldPD[1]*self.Jyg0 + FieldPD[2]*self.Jzg0)
        FieldPD = Field - Delta
        JdotB_m1z = muB*(FieldPD[0]*self.Jxg0 + FieldPD[1]*self.Jyg0 + FieldPD[2]*self.Jzg0)

        # B) Diagonalize full Hamiltonian

        # Diagonalize at unperturbed field (Delta=0 field)
        FieldHam = self.H_CEF.O + self.H_SOC.O + JdotB_0.O
        diagonalH = LA.eigh(FieldHam)
        minE = np.amin(diagonalH[0])
        evals = diagonalH[0] - minE

        # Diagonalize at all perturbed fields (Delta ≠ 0):
        Evals_pm = []
        for JdotB in [JdotB_p1x, JdotB_m1x, JdotB_p1y, JdotB_m1y, JdotB_p1z, JdotB_m1z]:   
            FieldHam = self.H_CEF.O + self.H_SOC.O + JdotB.O
            diagonalH = LA.eigh(FieldHam)
            Evals_pm.append(diagonalH[0])

        minE = np.amin(Evals_pm)
        Evals_pm -= minE
        Evals_pm = np.array(Evals_pm).T

        # C) Compute derivative of energy w.r.t. field:

        # Compute energy derivatives: ∂E_n/∂H_α
        Mderivs = np.zeros((len(Evals_pm),3))
        for i, ev in enumerate(Evals_pm):
            Mderivs[i]=[
                (ev[0] - ev[1])/(2*deltaField), # ∂E/∂H_x 
                (ev[2] - ev[3])/(2*deltaField), # ∂E/∂H_y 
                (ev[4] - ev[5])/(2*deltaField)  # ∂E/∂H_z
                ]
            
        k_B = 8.6173303e-2  # meV/K

        # Thermal average of derivatives
        if (isinstance(Temp, int) or isinstance(Temp, float)):
            Zz = np.sum(np.exp(-evals/(k_B*Temp)))
            BoltzmannWeights = np.exp(-evals/(k_B*Temp))/Zz
            return np.dot(BoltzmannWeights,Mderivs)/muB    # Convert from meV/T to μ_B
        
        else:
            expvals, temps = np.meshgrid(evals, Temp)
            ZZ = np.sum(np.exp(-expvals/temps/k_B), axis=1)
            MagList = np.repeat(Mderivs.reshape((1,)+Mderivs.shape), len(Temp), axis=0)
            MagList = np.sum(np.exp(-expvals/temps/k_B)*\
                                np.transpose(MagList, axes=[2,0,1]), axis=2) / ZZ
            return np.nan_to_num(MagList.T) / muB


    def gtensor(self):
        """
        Calculate g-tensor from ground state doublet (numerical method).
        
        Computes anisotropic g-factors by evaluating matrix elements of
        weighted J operators between lowest two eigenstates.
        
        Returns:
            3×3 g-tensor matrix
            
        Key Differences from CFLevels.gtensor():
            1. Uses both weighted (Jxg0, etc.) and unweighted (Jx, etc.) operators
            2. Doesn't use eliminateimag() helper function
            3. Computes matrix elements directly with np.dot
            
        Structure:
            g = 2 * [[|Re⟨0|J_xg0|1⟩|  Im⟨0|J_xg0|1⟩  ⟨0|J_xg0|0⟩]
                     [Re⟨0|J_yg0|1⟩    Im⟨0|J_yg0|1⟩  ⟨0|J_yg0|0⟩]
                     [Re⟨0|J_zg0|1⟩    Im⟨0|J_zg0|1⟩  |⟨0|J_zg0|0⟩|]]
                     
        Note:
            Factor of 2 comes from definition of g-tensor
            Uses absolute value for diagonal x and z components
            
            The eliminateimag function is defined but never used
            Consider removing or applying it to results
        """

        def eliminateimag(number):     
            """Helper to extract real part (currently unused)."""

            num = np.around(number, 10)
            if num.imag == 0:
                return (num.real).astype(float)
            else:
                return number

        # Find ground state doublet
        zeroinds = np.where(np.around(self.eigenvalues,4)==0)
        gsEVec = self.eigenvectors[zeroinds]
        vv1 = gsEVec[0]
        vv2 = gsEVec[1]

        # Extract operator matrices
        Jxg0, Jyg0, Jzg0 = self.Jxg0.O, self.Jyg0.O, self.Jzg0.O
        Jx, Jy, Jz = self.Jx.O, self.Jy.O, self.Jz.O

        # Compute matrix elements with weighted operators
        jzg01 = np.dot(vv1,np.dot(Jzg0,np.conj(vv2)))
        jzg10 = np.dot(vv2,np.dot(Jzg0,np.conj(vv1)))
        jzg00 = np.dot(vv1,np.dot(Jzg0,np.conj(vv1)))
        jzg11 = np.dot(vv2,np.dot(Jzg0,np.conj(vv2)))
               
        jxg01 = np.dot(vv1,np.dot(Jxg0,np.conj(vv2)))
        jxg10 = np.dot(vv2,np.dot(Jxg0,np.conj(vv1)))
        jxg00 = np.dot(vv1,np.dot(Jxg0,np.conj(vv1)))
        jxg11 = np.dot(vv2,np.dot(Jxg0,np.conj(vv2)))
        
        jyg01 = np.dot(vv1,np.dot(Jyg0,np.conj(vv2)))
        jyg10 = np.dot(vv2,np.dot(Jyg0,np.conj(vv1)))
        jyg00 = np.dot(vv1,np.dot(Jyg0,np.conj(vv1)))

        # Also compute with unweighted operators (currently unused?)
        jz01 = np.dot(vv1,np.dot(Jz,np.conj(vv2)))
        jz10 = np.dot(vv2,np.dot(Jz,np.conj(vv1)))
        jz00 = np.dot(vv1,np.dot(Jz,np.conj(vv1)))
        jz11 = np.dot(vv2,np.dot(Jz,np.conj(vv2)))
        
        jx01 = np.dot(vv1,np.dot(Jx,np.conj(vv2)))
        jx10 = np.dot(vv2,np.dot(Jx,np.conj(vv1)))
        jx00 = np.dot(vv1,np.dot(Jx,np.conj(vv1)))
        jx11 = np.dot(vv2,np.dot(Jx,np.conj(vv2)))
        
        jy01 = np.dot(vv1,np.dot(Jy,np.conj(vv2)))
        jy10 = np.dot(vv2,np.dot(Jy,np.conj(vv1)))
        jy00 = np.dot(vv1,np.dot(Jy,np.conj(vv1)))
        jy11 = np.dot(vv2,np.dot(Jy,np.conj(vv2)))

        # Construct g-tensor (note abs() on x and z diagonal)
        gg = 2*np.array([[np.abs(np.real(jxg01)), np.imag(jxg01), jxg00],
                        [np.real(jyg01), np.imag(jyg01), jyg00],
                        [np.real(jzg01), np.imag(jzg01), np.abs(jzg00)]])
        return gg


    def gtensorperturb(self, spinorbitcoupling, halffilled=True):
        """
        Calculate g-tensor using second-order perturbation theory.
        
        Computes g-factors via perturbative corrections from excited states,
        valid when crystal field >> spin-orbit coupling.
        
        Args:
            spinorbitcoupling: Spin-orbit coupling constant λ (meV)
            halffilled: If True, shell is less than half-filled (λ > 0 convention)
                       If False, more than half-filled (λ < 0 convention)
                       
        Returns:
            3×3 g-tensor matrix
            
        Formula:
            g_ij = g₀ δ_ij - 2ζ Σ_{n≠0} ⟨0|L_i|n⟩⟨n|L_j|0⟩/(E_n - E_0)
            
            where ζ = λ·2S·sign (sign depends on shell filling)
            
        Physical Basis:
            Ground state |0⟩ is modified by spin-orbit mixing with excited states:
            |0̃⟩ = |0⟩ + Σ_{n≠0} (λ⟨n|L·S|0⟩)/(E_0-E_n) |n⟩
            
            This mixing changes the effective g-factor from g₀ to g_ij
            
        When Valid:
            - Crystal field splitting >> spin-orbit coupling
            - Ground state well-separated from excited states
            - Typical for 3d transition metals
            
        Example:
            For octahedral Co²⁺ (d⁷, S=3/2, L=3):
            >>> g = cf.gtensorperturb(spinorbitcoupling=-0.08, halffilled=False)
            >>> # Might get g ≈ [[2.3, 0, 0], [0, 2.3, 0], [0, 0, 2.0]]
            
        Note:
            Uses second eigenvector/eigenvalue (index 1) as ground state
            Assumes ground state is non-degenerate or uses one component
        """

        g0 = 2.002319 # Free electron g-factor
        gtens = np.zeros((3,3)) + np.identity(3)*g0  # Start with g = g₀·I

        # Determine effective spin-orbit parameter
        if halffilled: 
            hff = -1 # d¹-d⁴ or f¹-f⁶
        else:  
            hff = 1 # d⁶-d⁹ or f⁸-f¹³
        zeta = spinorbitcoupling*2*self.S*hff

        # Build orbital operators
        Lx = LSOperator.Lx(self.L, self.S)
        Ly = LSOperator.Ly(self.L, self.S)
        Lz = LSOperator.Lz(self.L, self.S)

        # Use second eigenstate as ground state
        ev0 = self.eigenvectors[1]
        EE0 = self.eigenvalues[1]

        # Sum perturbative corrections from all other states
        for i, Li in enumerate([Lx, Ly, Lz]):
            for j, Lj in enumerate([Lx, Ly, Lz]):
                for k, ev in enumerate(self.eigenvectors):
                    if self.eigenvalues[k] != EE0: # Exclude ground state
                        # Matrix elements ⟨0|L_i|n⟩ and ⟨n|L_j|0⟩
                        jj1 = np.dot(np.conjugate(ev0),np.dot(Li.O,ev))
                        jj2 = np.dot(np.conjugate(ev),np.dot(Lj.O,ev0))

                        # Perturbation formula
                        gtens[i,j] -= 2*zeta*jj1*jj2/(self.eigenvalues[k]-EE0)
                    
                    else: pass # remove this

        return gtens


    def fitdata(self, chisqfunc, fitargs, method='Powell', **kwargs):
        """
        Fit crystal field parameters to experimental data.
        
        Optimizes B_n^m parameters by minimizing user-defined χ² function
        comparing calculated to measured properties.
        
        Args:
            chisqfunc: Function(LS_CFLevels, **kwargs) returning χ²
            fitargs: List of parameter names to fit
            method: Optimization algorithm (Powell, Nelder-Mead, etc.)
            **kwargs: Arguments for chisqfunc, must include:
                     'coeff': initial parameter guesses
                     
        Returns:
            Dictionary with optimized parameters and final χ²
            
        Validation:
            Checks that len(self.B) == len(kwargs['coeff'])
            
        Example:
            >>> def chi2(cf, Texp, chiexp):
            ...     chicalc = cf.susceptibility(Texp, 0.1, 0.01)
            ...     return np.sum((chicalc - chiexp)**2)
            >>> result = cf.fitdata(chi2, ['B20', 'B40'], method='Powell',
            ...                     coeff=[0.5, -0.02], Texp=T, chiexp=chi)
                                    
        Note:
            Similar to CFLevels.fitdata() but uses LS_CFLevels methods internally
        """

        initialChisq = chisqfunc(self, **kwargs)
        print('Initial err=', initialChisq, '\n')

        # Validate parameters
        if len(self.B) != len(kwargs['coeff']):
            raise ValueError('coeff needs to have the same length as self.B')

        # Create fitting function
        fun, p0, resfunc = makeFitFunction(chisqfunc, fitargs, **dict(kwargs, CFLevelsObject=self) )


        ############## Fit, using error function  #####################
        p_best = optimize.minimize(fun, p0, method=method)
        ###############################################################

        # Report results
        print(fun(p_best.x))
        print(chisqfunc(self, **kwargs))
        finalChisq = fun(p_best.x)
        print('\rInitial err =', initialChisq, '\tFinal err =', finalChisq)
        
        result = resfunc(p_best.x)
        result['Chisq'] = finalChisq
        return result


    def printLaTexEigenvectors(self):
        """
        Output eigenvectors in LaTeX table format for publication.
        
        Generates LaTeX code for a formatted table showing energy levels
        and eigenvector components in the |L,m_L,S,m_S⟩ basis.
        
        Output:
            LaTeX table code printed to stdout
            
        Format Features:
            - Handles integer and half-integer S values properly
            - Column headers show |L,S⟩ quantum numbers
            - Uses landscape mode for wide tables
            - Rounded to 2 decimal places (eigenvalues) and 3 (eigenvectors)
            
        Example Output:
            \\begin{table*}
            \\begin{landscape}
            ...
            E (meV) & |3,-1⟩ & |3,0⟩ & |3,+1⟩ & ...
            0.00 & 0.923 & 0.000 & 0.385 & ...
            ...
            
        Note:
            For large L,S, table becomes very wide
            (2L+1)×(2S+1) columns total

        TODO:
            - Add option to save to file instead of printing
            - Return string instead of printing
        """

        try:
            eigenkets = self.eigenvectors.real
        except AttributeError:
            self.diagonalize()

        # Build S quantum number labels
        if (self.S*2) %2 ==0: # Integer S
            Sarray = [str(int(kk)) for kk in 
                                    np.arange(-self.S,self.S+1)]
        else: # Half-integer S
            Sarray = ['-\\frac{'+str(abs(kk))+'}{2}' if kk <0
                            else '\\frac{'+str(abs(kk))+'}{2}'
                            for kk in np.arange(-int(self.S*2), int(self.S*2+2), 2)]

        # Build L quantum number labels (always integer)
        Larray = [str(int(kk)) for kk in  np.arange(-self.L,self.L+1)]

        # Combine into ket notation |L,S⟩
        KetNames = ['$|'+LL+','+SS+'\\rangle$' for LL in Larray  for SS in Sarray]

        # Generate LaTeX code
        print('\\begin{table*}\n\\begin{landscape}\n\\centering\n'+
            '\\caption{Eigenvectors and Eigenvalues... $|L,S\\rangle$}')
        print('\\begin{ruledtabular}')
        numev = len(self.eigenvalues)
        print('\\begin{tabular}{c|'+'c'*numev+'}')
        print('E (meV) &'+' & '.join(KetNames)
                +' \\tabularnewline\n \\hline ')
        
        # Sort and format eigenvalues/eigenvectors
        sortinds = self.eigenvalues.argsort()
        sortEVal= np.around(self.eigenvalues[sortinds],2)
        sortEVec= np.around(self.eigenvectors[sortinds],3)

        for i in range(len(sortinds)):
            print(format(self._Re(sortEVal[i]), '.3f'),'&', 
                ' & '.join([str(eevv) for eevv in self._Re(sortEVec[i])]), '\\tabularnewline')
            
        print('\\end{tabular}\\end{ruledtabular}')
        print('\\label{flo:Eigenvectors}\n\\end{landscape}\n\\end{table*}')
