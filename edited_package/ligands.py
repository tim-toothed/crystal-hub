"""
Point charge model calculations for crystal field theory.

Provides classes for computing crystal field parameters from ligand positions
using the point charge approximation. Supports both J-basis (rare earth) and
LS-coupling (transition metal) schemes.

Key Features:
- Automatic geometry handling from CIF files
- Point charge model with customizable charges
- Support for symmetry-equivalent ligands
- Parameter fitting to experimental data
- Coordinate system transformations

Classes:
    Ligands: Point charge calculations in J-basis (rare earth ions)
    LS_Ligands: Point charge calculations in LS-coupling basis (transition metals)
"""

import numpy as np
from scipy import optimize

from .lattice_class import lattice
from .constants import calculate_tesseral_harmonic, theta, calculate_radial_integral_RE, Constant, LStheta, PFalpha, PFbeta, calculate_radial_integral_TM, ION_NUMS_RARE_EARTH, is_half_filled
from .stevens_operators import StevensOp, LS_StevensOp
from .create_fit_function import makeFitFunction
from .cf_levels import CFLevels, LS_CFLevels
from .fundamental_operators import LSOperator


class Ligands:
    """
    Point charge model for rare earth ions in J-basis.
    
    Computes crystal field parameters from ligand positions using point charge
    approximation. Automatically handles coordinate transformations and Stevens
    operator construction.
    
    Attributes:
        ion (str): Ion symbol (e.g., 'Yb3+', 'Er3+')
        bonds (ndarray): Ligand position vectors from ion in Cartesian coordinates (Å)
        bondlen (ndarray): Bond lengths (Å)
        latt (lattice): Unit cell for coordinate transformations
        H (ndarray): Crystal field Hamiltonian matrix
        B (ndarray): Stevens parameter values
        IonCharge (int): Central ion charge
        suppressmm (bool): Whether to suppress negative m terms
    
    Example:
        >>> positions = np.array([[1,0,0], [-1,0,0], [0,1,0], [0,-1,0]])
        >>> lig = Ligands('Yb3+', positions)
        >>> cf = lig.PointChargeModel(LigandCharge=-2)
        >>> cf.diagonalize()
    """

    def __init__(self,ion,ligandPos, latticeParams=None, ionPos=[0,0,0]):
        """
        Initialize ligand geometry around magnetic ion.
        
        Args:
            ion: Ion symbol (e.g., 'Yb3+', 'Dy3+')
            ligandPos: Ligand positions as list/array of 3D vectors
                Can be fractional (ABC) or Cartesian depending on latticeParams
            latticeParams: [a,b,c,α,β,γ] unit cell parameters (Angstroms, degrees)
                If None, assumes Cartesian coordinates (cubic cell)
            ionPos: Central ion position in same coordinate system as ligandPos
        
        Raises:
            LookupError: If latticeParams doesn't have exactly 6 components
        """

        lp = latticeParams
        if lp == None:
            self.latt = lattice(1,1,1,90,90,90)
        elif len(lp) != 6:
            raise LookupError("latticeParams needs to have 6 components: a,b,c,alpha,beta,gamma")
        else:
            self.latt = lattice(lp[0], lp[1], lp[2], lp[3], lp[4], lp[5])

        self.bonds = np.array([O - np.array(ionPos) for O in ligandPos])
        self.bonds = self.latt.cartesian(self.bonds).astype('float')
        self.bondlen = np.linalg.norm(self.bonds, axis=1)
        self.ion = ion

    def rotateLigands(self, oldaxis, newaxis):
        """
        Rotate ligand geometry to align oldaxis with newaxis.
        
        Performs rotation so that a crystallographic direction (oldaxis)
        becomes aligned with a new reference direction (newaxis). Useful
        for standardizing orientation before calculations.
        
        Args:
            oldaxis: Original direction vector (e.g., [1,1,0])
            newaxis: Target direction vector (e.g., [0,0,1])
        
        Note:
            Modifies self.bonds in-place using Rodrigues rotation formula.
        """
        rotationAxis = np.cross(newaxis,oldaxis)
        rotationAngle = np.arccos(np.dot(newaxis,oldaxis)/(np.linalg.norm(newaxis)*np.linalg.norm(oldaxis)))
        self.bonds = np.array([self._rotateMatrix(b,rotationAxis,rotationAngle) for b in self.bonds])

    def rotateLigandsZ(self, oldaxis):
        """
        Rotate ligands around z-axis to place oldaxis along x-axis.
        
        Secondary rotation after rotateLigands() to fix azimuthal orientation.
        Places specified axis in xy-plane along x direction.
        
        Args:
            oldaxis: Direction vector to align with x-axis (z-component ignored)
        """
        zrotation = np.arctan(oldaxis[1]/oldaxis[0])
        self.bonds = np.array([self._rotateMatrix(b,np.array([0,0,1]),-zrotation) for b in self.bonds])


    def _rotateMatrix(self,matrixin,axis,angle):
        """
        Rotate vector/matrix using Rodrigues rotation formula.
        
        Internal method for arbitrary axis rotation. Handles both vectors
        (simple rotation) and matrices (similarity transform R·M·R^T).
        
        Args:
            matrixin: 3-element vector or 3×3 matrix to rotate
            axis: Rotation axis vector (doesn't need to be normalized)
            angle: Rotation angle in radians
        
        Returns:
            Rotated vector (3-element) or matrix (3×3)
        
        Reference:
            http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
        """
        u, v, w = axis[0], axis[1], axis[2]
        norm = u**2 + v**2 + w**2
        
        rotmatrix = np.zeros((3,3))
        rotmatrix[0,0] = (u**2 +(v**2 + w**2)*np.cos(angle)) / norm
        rotmatrix[0,1] = (u*v*(1- np.cos(angle)) - w*np.sqrt(norm)*np.sin(angle)) / norm
        rotmatrix[0,2] = (u*w*(1- np.cos(angle)) + v*np.sqrt(norm)*np.sin(angle)) / norm
        rotmatrix[1,0] = (u*v*(1- np.cos(angle)) + w*np.sqrt(norm)*np.sin(angle)) / norm
        rotmatrix[1,1] = (v**2 +(u**2 + w**2)*np.cos(angle)) / norm
        rotmatrix[1,2] = (v*w*(1- np.cos(angle)) - u*np.sqrt(norm)*np.sin(angle)) / norm
        rotmatrix[2,0] = (u*w*(1- np.cos(angle)) - v*np.sqrt(norm)*np.sin(angle)) / norm
        rotmatrix[2,1] = (v*w*(1- np.cos(angle)) + u*np.sqrt(norm)*np.sin(angle)) / norm
        rotmatrix[2,2] = (w**2 +(v**2 + u**2)*np.cos(angle)) / norm

        # Simple matrix multiplication of matrixin is a vector
        if matrixin.size == 3:
            return np.dot(rotmatrix, matrixin)
        # R*m*R^T if matrixin is a matrix
        elif matrixin.size == 9:
            return np.dot(rotmatrix, np.dot(matrixin, rotmatrix.transpose() ))

    # <- повторяет exportLigandCif
    def exportCif(self, filename):
        """
        See exportLigandCif
        """
        exportLigandCif(self, filename)


    def PointChargeModel(self, symequiv=None, LigandCharge=-2,IonCharge=3, printB = True, 
                            suppressminusm = False, ionL=None):
        """
        Calculate crystal field Hamiltonian via point charge model.
        
        Computes Stevens operator coefficients B_n^m by treating ligands as
        point charges. Uses tesseral harmonics and radial integrals from
        literature values (Hutchings 1964, Freeman & Desclaux 1979).
        
        Args:
            symequiv: Symmetry equivalence list for grouped ligand charges
                If None, all ligands use same charge (or list of charges)
            LigandCharge: Ligand charge(s) in units of |e|
                - float: Same charge for all ligands
                - list: Individual charges per ligand (length must match bonds)
                - list with symequiv: Charges indexed by symmetry group
            IonCharge: Central ion charge in units of |e| (default 3 for RE3+)
            printB: If True, print calculated B_n^m parameters to console
            suppressminusm: If True, only include m≥0 terms (reduces Hamiltonian size)
            ionL: Override angular momentum (default: use database value for ion)
        
        Returns:
            CFLevels object with Hamiltonian set, containing:
                - H: Full Hamiltonian matrix
                - B: List of non-zero Stevens coefficients (meV)
                - O: List of corresponding Stevens operator matrices
                - BnmLabels: Parameter labels like 'B_4^2'
        
        Raises:
            ValueError: If LigandCharge list length doesn't match number of bonds
        
        Example:
            >>> lig = Ligands('Yb3+', ligand_positions)
            >>> # Octahedral: 6 O2- ligands with same charge
            >>> cf = lig.PointChargeModel(LigandCharge=-2)
            >>> # Lower symmetry: different charges
            >>> cf = lig.PointChargeModel(LigandCharge=[-2,-2,-1.5,-1.5])
        
        Note:
            - Energy scale set by ahc = α·ℏ·c = 1.43996e4 meV·Å
            - Bohr radius a₀ = 0.52917721067 Å
            - Includes all n=2,4,6 Stevens operators up to m=±n
        """

        self.IonCharge = IonCharge
        # Lock suppressmm into whatever it was when PointChargeModel was first called.
        try: self.suppressmm
        except AttributeError:
            self.suppressmm = suppressminusm


        if symequiv == None:
            try:
                if len(LigandCharge) == len(self.bonds):
                    charge = LigandCharge
                else:
                    charge = [LigandCharge]*len(self.bonds)
            except TypeError:
                charge = [LigandCharge]*len(self.bonds)

        else:
            charge = [0]*len(self.bonds)
            for i,se in enumerate(symequiv):
                charge[i] = LigandCharge[se]

        
        ion=self.ion
        if ionL == None:
            ionJ = ION_NUMS_RARE_EARTH[ion][2]
        else: ionJ = ionL

        ahc = 1.43996e4  #Constant to get the energy in units of meV = alpha*hbar*c
        a0 = 0.52917721067    #Bohr radius in \AA

        self.H = np.zeros((int(2*ionJ+1), int(2*ionJ+1)),dtype = complex)
        self.B = []
        OOO = []
        nonzeroB = []
        bnm_labels = []

        if self.suppressmm == False:  nmrange = [[n,m] for n in range(2,8,2) for m in range(-n,n+1)]
        elif self.suppressmm == True:   nmrange = [[n,m] for n in range(2,8,2) for m in range(0,n+1)]

        for n,m in nmrange:
            # 1)  Compute gamma
            gamma = 0
            for i in range(len(self.bonds)):
                #print(np.squeeze(charge[i]))
                gamma += 4*np.pi/(2*n+1)*np.squeeze(charge[i]) *\
                            calculate_tesseral_harmonic(n,m, self.bonds[i][0], self.bonds[i][1], self.bonds[i][2])/\
                            (self.bondlen[i]**(n+1))

            # 2)  Compute CEF parameter
            B = -gamma * ahc* a0**n * Constant(n,m) * calculate_radial_integral_RE(ion,n) * theta(ion,n)
            if printB ==True: print('B_'+str(n),m,' = ',np.around(B,decimals=8))
            if np.around(B,decimals=7) != 0:
                OOO.append(StevensOp(ionJ,n,m))
                nonzeroB.append(B)
                bnm_labels.append('B_{}^{}'.format(n,m))

            if np.around(B,decimals=9) != 0:
                self.H += B*StevensOp(ionJ,n,m)
            self.B.append(B)
        self.B = np.array(self.B)

        newobj = CFLevels.Hamiltonian(self.H)
        newobj.O = OOO
        newobj.B = nonzeroB
        newobj.BnmLabels = bnm_labels
        newobj.ion = self.ion
        return newobj

    # <- Legacy - убрать
    def FitChargesNeutrons(self, chisqfunc, fitargs, method='Powell', **kwargs):
        """
        Legacy alias for FitCharges(). Kept for backward compatibility.
        
        Args/Returns: See FitCharges()
        """
        return self.FitCharges(chisqfunc, fitargs, method='Powell', **kwargs)

    def FitCharges(self, chisqfunc, fitargs, method='Powell', **kwargs):
        """
        Fit ligand charges to experimental data.
        
        Optimizes ligand charges to minimize chi-squared against experimental
        observables (susceptibility, magnetization, neutron intensities, etc.).
        
        Args:
            chisqfunc: Chi-squared function with signature:
                chisqfunc(CFLevels_object, **kwargs) -> float
            fitargs: Dictionary of initial parameter guesses, e.g.:
                {'LigandCharge': [-2.0, -1.8]}
            method: scipy.optimize.minimize method ('Powell', 'Nelder-Mead', etc.)
            **kwargs: Additional arguments passed to chi-squared function
                Common: symequiv, temperature arrays, field arrays
        
        Returns:
            tuple: (newH, finalvals)
                - newH: CFLevels object with optimized Hamiltonian
                - finalvals: Dictionary of optimized parameters
        
        Example:
            >>> def chi2(cf, temps, data):
            ...     chi_calc = cf.susceptibility(temps, Field=0.1)
            ...     return np.sum((chi_calc - data)**2)
            >>> 
            >>> initial = {'LigandCharge': [-2.0]}
            >>> cf_fit, params = lig.FitCharges(chi2, initial, 
            ...                                  temps=temps, data=chi_exp)
        
        Note:
            Prints fitting progress, initial/final chi-squared, and optimized
            Stevens parameters to console.
        """

        # Define function to be fit
        fun, p0, resfunc = makeFitFunction(chisqfunc, fitargs, **dict(kwargs, LigandsObject=self) )

        print('\tFitting...')
        ############## Fit, using error function  #####################
        p_best = optimize.minimize(fun, p0, method=method)
        ###############################################################

        try:
            initialChisq, finalChisq = fun(p0), fun(p_best.x)           # <- неиспользуемые переменные
            finalvals = resfunc(p_best.x)
        except IndexError:
            initialChisq, finalChisq = fun(p0), fun([float(p_best.x)])  # <- неиспользуемые переменные
            finalvals = resfunc([float(p_best.x)])

        # split back into values
        finalCharges = finalvals['LigandCharge']

        # Print results
        print("\n#*********************************")
        print("# Final Stevens Operator Values")
        try:
            newH = self.PointChargeModel(kwargs['symequiv'], finalCharges, printB=True)
        except KeyError:
            print(float(p_best.x))
            newH = self.PointChargeModel(LigandCharge=[float(p_best.x)], printB=True)
        newH.diagonalize()
        print("\nFinal Charges: ", finalCharges)
        print('Final EigenValues: ', np.around(np.sort(newH.eigenvalues.real),3))

        return newH, finalvals

def exportLigandCif(ligands, filename):
    """
    Export Ligands object to CIF file showing local environment.

    Utility function to write minimal CIF containing only the central ion
    and its nearest-neighbor ligands for visualization.

    Args:
        ligands: Ligands object instance
        filename: Output path (.cif extension added if missing)

    Note:
        Creates P1 space group with 10Å cubic cell centered at (0.5,0.5,0.5).
        Ligands labeled as 'S2-' regardless of actual chemistry.
    """

    if not filename.endswith('.cif'):
        filename = filename + '.cif' 

    with open(filename, 'w') as f:
        f.write('# Output from PyCrystalField showing the ligand environment\n\n')
        f.write('loop_\n'+\
				'_publ_author_name\n'+\
				"'Someone, A.'\n"+\
				"'Someone, B.'\n"+
				'_cell_length_a 10.0\n'+\
				'_cell_length_b 10.0\n'+\
				'_cell_length_c 10.0\n'+\
				'_cell_angle_alpha 90.\n'+\
				'_cell_angle_beta 90.\n'+\
				'_cell_angle_gamma 90.\n'+\
				'_cell_volume 1000.0\n'+\
				"_symmetry_space_group_name_H-M 'P 1'\n"+\
				'loop_\n'+\
				'_symmetry_equiv_pos_site_id\n'+\
				'_symmetry_equiv_pos_as_xyz\n'+\
				"1 'x, y, z'\n"+\
				'loop_\n'+\
				'_atom_type_symbol\n'+\
				ligands.ion + '\n'+\
				'S2-\n'+\
				'loop_\n'+\
				'_atom_site_label\n'+\
				'_atom_site_fract_x\n'+\
				'_atom_site_fract_y\n'+\
				'_atom_site_fract_z\n'+\
				'_atom_site_occupancy\n'+\
				ligands.ion +' 0.5 0.5 0.5 1. \n')
        for b in ligands.bonds:
            f.write('S1 '+ ' '.join([str(bi/10+0.5) for bi in b])+ ' 1. \n')

class LS_Ligands:
    """
    Point charge model for transition metals with LS-coupling basis.

    LS basis functions are a set of quantum mechanical wave functions
    that describe the electronic states of an atom 
    in the L-S (Russell-Saunders) coupling scheme. 
    This model is used for lighter atoms where 
    the electrostatic repulsion between electrons
    is stronger than the spin-orbit interaction.

    For 3d transition metal ions where spin-orbit coupling is weak compared
    to crystal field splitting. Builds Hamiltonian in |L,m_L⟩⊗|S,m_S⟩ basis
    and adds spin-orbit coupling λ·L·S as perturbation.
    
    Attributes:
        ion (str): Ion symbol or custom label
        ionL (float): Orbital angular momentum quantum number
        ionS (float): Spin quantum number
        bonds (np.ndarray): Ligand vectors in Cartesian coordinates (Å)
        bondlen (np.ndarray): Bond lengths (Å)
        H_SOC (LSOperator): Spin-orbit coupling Hamiltonian
        H_CEF (LSOperator): Crystal field Hamiltonian (after PointChargeModel)
        B (np.ndarray): Stevens coefficients (meV)
    
    Example:
        >>> # Ni2+: 3d8, L=3, S=1
        >>> lig = LS_Ligands('Ni2+', positions, SpinOrbitCoupling=-315)
        >>> cf = lig.PointChargeModel(LigandCharge=-2)
        >>> cf.diagonalize()
    """

    def __init__(self,ion, ligandPos, SpinOrbitCoupling, latticeParams=None, ionPos=[0,0,0]):
        """
        Initialize LS-basis ligand calculation.
        
        Args:
            ion: Either ion symbol string (e.g., 'Ni3+') or list [label, S, L]
                If list: ['MyIon', 1.0, 2.0] specifies custom L and S
            ligandPos: Ligand position vectors (fractional or Cartesian)
            SpinOrbitCoupling: λ parameter in meV for H_SOC = λ·L·S
                Typically negative for less-than-half-filled shells
            latticeParams: [a,b,c,α,β,γ] unit cell (None → Cartesian)
            ionPos: Central ion position
        
        Note:
            For standard ions, L and S are retrieved from ION_NUMS_RARE_EARTH
            database. For custom ions, provide [label, S, L] list.
        """
        lp = latticeParams
        if lp == None:
            self.latt = lattice(1,1,1,90,90,90)
        elif len(lp) != 6:
            raise LookupError("latticeParams needs to have 6 components: a,b,c,alpha,beta,gamma")
        else:
            self.latt = lattice(lp[0], lp[1], lp[2], lp[3], lp[4], lp[5])

        self.bonds = np.array([np.array(O) - np.array(ionPos) for O in ligandPos])
        self.bonds = self.latt.cartesian(self.bonds).astype('float')
        self.bondlen = np.linalg.norm(self.bonds, axis=1)

        if isinstance(ion, str):
            self.ion = ion
            self.ionS = ION_NUMS_RARE_EARTH[ion][0]
            self.ionL = ION_NUMS_RARE_EARTH[ion][1]
        else:
            self.ion = ion[0]
            self.ionS = ion[1]
            self.ionL = ion[2]

        # Now, define the spin orbit coupling (so we don't have to re-define it 
        # every time we build the point charge model).
        Sx = LSOperator.Sx(self.ionL, self.ionS)
        Sy = LSOperator.Sy(self.ionL, self.ionS)
        Sz = LSOperator.Sz(self.ionL, self.ionS)
        Lx = LSOperator.Lx(self.ionL, self.ionS)
        Ly = LSOperator.Ly(self.ionL, self.ionS)
        Lz = LSOperator.Lz(self.ionL, self.ionS)

        self.H_SOC = Lx*Sx + Ly*Sy + Lz*Sz
        LdotS = self.H_SOC.O*1.0
        if np.sum(LdotS.imag) == 0: LdotS = LdotS.real
        self.H_SOC.O = SpinOrbitCoupling*LdotS

    def rotateLigands(self, oldaxis, newaxis):
        """Rotate geometry to align oldaxis with newaxis. See Ligands.rotateLigands()."""
        rotationAxis = np.cross(newaxis,oldaxis)
        rotationAngle = np.arccos(np.dot(newaxis,oldaxis)/(np.linalg.norm(newaxis)*np.linalg.norm(oldaxis)))
        self.bonds = np.array([self._rotateMatrix(b,rotationAxis,rotationAngle) for b in self.bonds])

    def rotateLigandsZ(self, oldaxis):
        """Rotate around z to place oldaxis along x. See Ligands.rotateLigandsZ()."""
        zrotation = np.arctan(oldaxis[1]/oldaxis[0])
        self.bonds = np.array([self._rotateMatrix(b,np.array([0,0,1]),-zrotation) for b in self.bonds])


    def _rotateMatrix(self,matrixin,axis,angle):
        """Rodrigues rotation. See Ligands._rotateMatrix()."""
        u, v, w = axis[0], axis[1], axis[2]
        norm = u**2 + v**2 + w**2
        
        rotmatrix = np.zeros((3,3))
        rotmatrix[0,0] = (u**2 +(v**2 + w**2)*np.cos(angle)) / norm
        rotmatrix[0,1] = (u*v*(1- np.cos(angle)) - w*np.sqrt(norm)*np.sin(angle)) / norm
        rotmatrix[0,2] = (u*w*(1- np.cos(angle)) + v*np.sqrt(norm)*np.sin(angle)) / norm
        rotmatrix[1,0] = (u*v*(1- np.cos(angle)) + w*np.sqrt(norm)*np.sin(angle)) / norm
        rotmatrix[1,1] = (v**2 +(u**2 + w**2)*np.cos(angle)) / norm
        rotmatrix[1,2] = (v*w*(1- np.cos(angle)) - u*np.sqrt(norm)*np.sin(angle)) / norm
        rotmatrix[2,0] = (u*w*(1- np.cos(angle)) - v*np.sqrt(norm)*np.sin(angle)) / norm
        rotmatrix[2,1] = (v*w*(1- np.cos(angle)) + u*np.sqrt(norm)*np.sin(angle)) / norm
        rotmatrix[2,2] = (w**2 +(v**2 + u**2)*np.cos(angle)) / norm

        # Simple matrix multiplication of matrixin is a vector
        if matrixin.size == 3:
            return np.dot(rotmatrix, matrixin)
        # R*m*R^T if matrixin is a matrix
        elif matrixin.size == 9:
            return np.dot(rotmatrix, np.dot(matrixin, rotmatrix.transpose() ))


    def PointChargeModel(self,  symequiv=None, LigandCharge=-2,IonCharge=1,
                        printB = True, suppressminusm = False):
        """
        Calculate LS-basis crystal field Hamiltonian.
        
        Computes crystal field in orbital basis, then expands to full LS
        space via tensor product with spin identity: H_CF ⊗ I_spin.
        
        Args:
            symequiv: Ligand charge grouping (None for uniform)
            LigandCharge: Charge(s) in units of |e|
            IonCharge: Central ion charge (default 1 for TM2+/3+)
            printB: Print Stevens coefficients
            suppressminusm: Use only m≥0 terms
        
        Returns:
            LS_CFLevels object with:
                - H = H_CEF + H_SOC combined Hamiltonian
                - Separate H_CEF and H_SOC accessible
        
        Note:
            Uses LStheta() for Stevens factors (different from J-basis theta).
            Only n=2,4,6 operators included (n=6 rare for 3d ions).
        """

        # Lock suppressmm into whatever it was when PointChargeModel was first called.
        self.IonCharge = IonCharge
        try: self.suppressmm
        except AttributeError:
            self.suppressmm = suppressminusm

        if symequiv == None:
            try:
                if len(LigandCharge) == len(self.bonds):
                    charge = LigandCharge
                else:
                    charge = [LigandCharge]*len(self.bonds)
            except TypeError:
                charge = [LigandCharge]*len(self.bonds)

        else:
            charge = [0]*len(self.bonds)
            for i,se in enumerate(symequiv):
                charge[i] = LigandCharge[se]
        
        ion=self.ion

        ahc = 1.43996e4  #Constant to get the energy in units of meV = alpha*hbar*c
        a0 = 0.52917721067    #Bohr radius in \AA

        H = np.zeros((int(2*self.ionL+1), int(2*self.ionL+1)),dtype = complex)
        self.B = []
        OOO = []
        nonzeroB = []

        self.H_nocharge = [[]]
        if self.suppressmm == False:  nmrange = [[n,m] for n in range(2,8,2) for m in range(-n,n+1)]
        elif self.suppressmm == True:   nmrange = [[n,m] for n in range(2,8,2) for m in range(0,n+1)]
        for n,m in nmrange:
            # 1)  Compute gamma
            gamma = 0
            for i in range(len(self.bonds)):

                gamma += 4*np.pi/(2*n+1)*charge[i] *\
                            calculate_tesseral_harmonic(n,m, self.bonds[i][0], self.bonds[i][1], self.bonds[i][2])/\
                            (self.bondlen[i]**(n+1))

            # 2)  Compute CEF parameter
            B = -gamma * ahc* a0**n * Constant(n,m) * calculate_radial_integral_RE(ion,n) * LStheta(ion,n)
            if printB ==True: print('B_'+str(n),m,' = ',np.around(B,decimals=8))
            if np.around(B,decimals=8) != 0:
                OOO.append(LS_StevensOp(self.ionL,self.ionS,n,m))
                nonzeroB.append(B)
            if np.around(B,decimals=10) != 0:
                H += B*StevensOp(self.ionL,n,m)
            self.B.append(B)
        self.B = np.array(self.B)

        # Convert Hamiltonian to full LS basis
        H_CEF_O = np.hstack(np.hstack(np.multiply.outer(H, np.identity(int(2*self.ionS+1)))))
        self.H_CEF = LSOperator(self.ionL, self.ionS)
        self.H_CEF.O = H_CEF_O

        newobj = LS_CFLevels.Hamiltonian(self.H_CEF, self.H_SOC, self.ionL, self.ionS)
        newobj.O = OOO
        newobj.B = nonzeroB
        return newobj


    def TMPointChargeModel(self, l=2, symequiv=None, LigandCharge= -2, IonCharge=1,
                        printB = True, suppressminusm = False):
        """
        Specialized point charge model for 3d transition metals.
        
        Uses Pott's factor corrections (PFalpha, PFbeta) and experimental
        radial integrals from literature. Automatically handles half-filled
        vs non-half-filled shell differences.
        
        Args:
            l: Orbital type (2 for d-orbitals, 3 for f-orbitals)
            symequiv: Charge grouping
            LigandCharge: Ligand charge(s)
            IonCharge: Central ion charge
            printB: Print parameters
            suppressminusm: m≥0 only
        
        Returns:
            LS_CFLevels object with TM-specific radial integrals
        
        Note:
            - Uses calculate_radial_integral_TM() from constants
            - Half-filled detection via is_half_filled(ion)
            - Limited to n=2,4 (n=6 negligible for 3d)
        """
        halffilled = is_half_filled(self.ion)
        
        self.IonCharge = IonCharge
        # Lock suppressmm into whatever it was when PointChargeModel was first called.
        try: self.suppressmm
        except AttributeError:
            self.suppressmm = suppressminusm


        if symequiv == None:
            try:
                if len(LigandCharge) == len(self.bonds):
                    charge = LigandCharge
                else:
                    charge = [LigandCharge]*len(self.bonds)
            except TypeError:
                charge = [LigandCharge]*len(self.bonds)

        else:
            charge = [0]*len(self.bonds)
            for i,se in enumerate(symequiv):
                charge[i] = LigandCharge[se]



        ahc = 1.43996e4  #Constant to get the energy in units of meV = alpha*hbar*c
        a0 = 0.52917721067    #Bohr radius in \AA

        H = np.zeros((int(2*self.ionL+1), int(2*self.ionL+1)),dtype = complex)
        self.B = []
        OOO = []
        nonzeroB = []

        TM_LStheta = {2: PFalpha(self.ionL,self.ionS,l,halffilled), 
                    4: PFbeta(self.ionL,self.ionS,l,halffilled)}

        self.H_nocharge = [[]]
        if self.suppressmm == False:  nmrange = [[n,m] for n in range(2,6,2) for m in range(-n,n+1)]
        elif self.suppressmm == True:   nmrange = [[n,m] for n in range(2,6,2) for m in range(0,n+1)]
        for n,m in nmrange:
            # 1)  Compute gamma
            gamma = 0
            for i in range(len(self.bonds)):

                gamma += 4*np.pi/(2*n+1)*charge[i] *\
                            calculate_tesseral_harmonic(n,m, self.bonds[i][0], self.bonds[i][1], self.bonds[i][2])/\
                            (self.bondlen[i]**(n+1))

            # 2)  Compute CEF parameter
            B = -gamma * ahc* a0**n * Constant(n,m) * calculate_radial_integral_TM(self.ion, n) * TM_LStheta[n]
            if printB ==True: print('B_'+str(n),m,' = ',np.around(B,decimals=8))
            if np.around(B,decimals=8) != 0:
                OOO.append(LS_StevensOp(self.ionL,self.ionS,n,m))
                nonzeroB.append(B)
            if np.around(B,decimals=10) != 0:
                H += B*StevensOp(self.ionL,n,m)
            self.B.append(B)
        self.B = np.array(self.B)

        # Convert Hamiltonian to full LS basis
        H_CEF_O = np.hstack(np.hstack(np.multiply.outer(H, np.identity(int(2*self.ionS+1)))))
        self.H_CEF = LSOperator(self.ionL, self.ionS)
        self.H_CEF.O = H_CEF_O

        #self.H = self.H_CEF + self.H_LS
        newobj = LS_CFLevels.Hamiltonian(self.H_CEF, self.H_SOC, self.ionL, self.ionS)
        newobj.O = OOO
        newobj.B = nonzeroB
        return newobj



    def UnknownTMPointChargeModel(self, radialintegrals, halffilled=True, l=2,
                        symequiv=None, LigandCharge= -2, IonCharge=1,
                        printB = True, suppressminusm = False):
        """
        Point charge model for TM ions not in standard database.
        
        For exotic oxidation states or 4d/5d ions where radial integrals
        aren't tabulated. User provides custom radial integral values.
        
        Args:
            radialintegrals: Dict mapping n→⟨r^n⟩ values, e.g., {2: 1.5, 4: 0.8}
            halffilled: True for d5/f7, False otherwise (affects Pott's factors)
            l: Orbital quantum number (2 for d, 3 for f)
            symequiv: Charge grouping
            LigandCharge: Ligand charge(s)
            IonCharge: Central ion charge
            printB: Print parameters
            suppressminusm: m≥0 only
        
        Returns:
            LS_CFLevels object using custom radial integrals
        
        Example:
            >>> # For hypothetical Tc4+ (4d3)
            >>> rad_int = {2: 2.1, 4: 1.3}  # Custom values
            >>> lig = LS_Ligands(['Tc4+', 1.5, 2], positions, lambda_soc=-200)
            >>> cf = lig.UnknownTMPointChargeModel(rad_int, halffilled=False)
        """

        self.IonCharge = IonCharge
        # Lock suppressmm into whatever it was when PointChargeModel was first called.
        try: self.suppressmm
        except AttributeError:
            self.suppressmm = suppressminusm

        if symequiv == None:
            try:
                if len(LigandCharge) == len(self.bonds):
                    charge = LigandCharge
                else:
                    charge = [LigandCharge]*len(self.bonds)
            except TypeError:
                charge = [LigandCharge]*len(self.bonds)

        else:
            charge = [0]*len(self.bonds)
            for i,se in enumerate(symequiv):
                charge[i] = LigandCharge[se]


        ahc = 1.43996e4  #Constant to get the energy in units of meV = alpha*hbar*c
        a0 = 0.52917721067    #Bohr radius in \AA

        H = np.zeros((int(2*self.ionL+1), int(2*self.ionL+1)),dtype = complex)
        self.B = []
        OOO = []
        nonzeroB = []

        TM_LStheta = {2: PFalpha(self.ionL,self.ionS,l,halffilled), 
                    4: PFbeta(self.ionL,self.ionS,l,halffilled)}

        self.H_nocharge = [[]]
        if self.suppressmm == False:  nmrange = [[n,m] for n in range(2,6,2) for m in range(-n,n+1)]
        elif self.suppressmm == True:   nmrange = [[n,m] for n in range(2,6,2) for m in range(0,n+1)]
        for n,m in nmrange:
            # 1)  Compute gamma
            gamma = 0
            for i in range(len(self.bonds)):

                gamma += 4*np.pi/(2*n+1)*charge[i] *\
                            calculate_tesseral_harmonic(n,m, self.bonds[i][0], self.bonds[i][1], self.bonds[i][2])/\
                            (self.bondlen[i]**(n+1))

            # 2)  Compute CEF parameter
            B = -gamma * ahc* a0**n * Constant(n,m) * radialintegrals[n] * TM_LStheta[n]
            if printB ==True: print('B_'+str(n),m,' = ',np.around(B,decimals=8))
            if np.around(B,decimals=8) != 0:
                OOO.append(LS_StevensOp(self.ionL,self.ionS,n,m))
                nonzeroB.append(B)
            if np.around(B,decimals=10) != 0:
                H += B*StevensOp(self.ionL,n,m)
            self.B.append(B)
        self.B = np.array(self.B)

        # Convert Hamiltonian to full LS basis
        H_CEF_O = np.hstack(np.hstack(np.multiply.outer(H, np.identity(int(2*self.ionS+1)))))
        self.H_CEF = LSOperator(self.ionL, self.ionS)
        self.H_CEF.O = H_CEF_O

        newobj = LS_CFLevels.Hamiltonian(self.H_CEF, self.H_SOC, self.ionL, self.ionS)
        newobj.O = OOO
        newobj.B = nonzeroB
        return newobj


    def FitChargesNeutrons(self, chisqfunc, fitargs, method='Powell', **kwargs):
        """
        Fit ligand charges to neutron scattering data.
        
        Optimizes charges to match experimental neutron intensities and
        transition energies. Similar to Ligands.FitCharges() but for
        LS-basis calculations.
        
        Args:
            chisqfunc: Chi-squared function taking LS_CFLevels object
            fitargs: Initial parameter dictionary
            method: Optimization method
            **kwargs: Additional args (symequiv, etc.)
        
        Returns:
            tuple: (optimized_LS_CFLevels, final_parameters)
        
        Note:
            Prints fitting results including eigenvalues and charges.
        """

        # Define function to be fit
        fun, p0, resfunc = makeFitFunction(chisqfunc, fitargs, **dict(kwargs, LigandsObject=self) )

        print('\tFitting...')
        ############## Fit, using error function  #####################
        p_best = optimize.minimize(fun, p0, method=method)
        ###############################################################

        initialChisq, finalChisq = fun(p0), fun(p_best.x) # <- неиспользуемые переменные

        # split back into values
        finalvals = resfunc(p_best.x)
        finalCharges = finalvals['LigandCharge']

        # Print results
        print("\n#*********************************")
        print("# Final Stevens Operator Values")
        newH = self.PointChargeModel(kwargs['symequiv'], finalCharges, printB=True)
        newH.diagonalize()
        print("\nFinal Charges: ", finalCharges)
        print('Final EigenValues: ', np.around(np.sort(newH.eigenvalues.real),3))

        return newH, finalvals
