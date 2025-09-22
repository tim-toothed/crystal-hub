import numpy as np
from optical_transitions import OpticalTransition
from pcf_lib.Operators import Ket, Operator
from constants import Jion
import numpy as np
from scipy import optimize
import scipy.linalg as LA
from scipy.special import wofz
from pcf_lib.form_factors import RE_FormFactor
from pcf_lib.CreateFitFunction import makeFitFunction
from pcf_lib.Operators import Ket, Operator
from pcf_lib.StevensOperators import StevensOp
from lande_g_factor import LandeGFactor

class CFLevels:
    """For calculating and fitting crystal field levels for an ion"""
    
    def __init__(self, StevensOperators, Parameters):
        """Initialize crystal field Hamiltonian and optical transition calculator"""
        # Create crystal field Hamiltonian
        self.H = np.sum([a*b for a, b in zip(StevensOperators, Parameters)], axis=0)
        self.O = StevensOperators  # Save Stevens operators for fitting
        self.B = Parameters        # Save crystal field parameters
        
        try:
            # Calculate total angular momentum quantum number
            self.J = (len(self.H) - 1.0) / 2.0
            
            # Initialize optical transition calculator
            self.opttran = OpticalTransition(
                Operator.Jx(self.J).O,        # J_x matrix
                Operator.Jy(self.J).O.imag,   # J_y matrix (real part)
                Operator.Jz(self.J).O         # J_z matrix
            )
        except (TypeError, AttributeError):
            # Handle cases where J cannot be determined
            pass

    @classmethod
    def Bdict(cls, ion, Bdict):
        ionJ = Jion[ion][-1]
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
        newcls = cls([0,0],[0,0])  # Create empty class so we can just define Hamiltonian
        newcls.H = Hamil
        newcls.J = (len(Hamil) -1.)/2
        newcls.opttran = OpticalTransition(Operator.Jx(newcls.J).O.real, Operator.Jy(newcls.J).O.imag, Operator.Jz(newcls.J).O.real)
        return newcls


    def newCoeff(self, newcoeff):
        self.B = np.array(newcoeff)
        newH = np.sum([a*b for a,b in zip(self.O, newcoeff)], axis=0)
        self.diagonalize(newH)

    def diagonalize(self, Hamiltonian=None, old=False):
        """A Hamiltonian can be passed to the function (used for data fits)
        or the initially defined hamiltonian is used."""
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
        '''same as above, but using the Scipy eig_banded function'''
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
        '''used in the diagonalize_banded function'''
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
        # for population factor weights
        beta = 1/(8.61733e-2*Temp)  # Boltzmann constant is in meV/K
        Z = sum([np.exp(-beta*en) for en in self.eigenvalues])
        # compute population factor
        pn = np.exp(-beta *self.eigenvalues[ii])/Z
        
        # compute amplitude
        mJn = self.opttran.transition(self.eigenvectors.real[ii] ,  self.eigenvectors.real[jj])
        return pn*mJn


    def neutronSpectrum(self, Earray, Temp, Ei, ResFunc, gamma = 0):
        # make angular momentum ket object
        #eigenkets = [Ket(ei) for ei in self.eigenvectors]

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
        '''calculate neutron spectrum with a custom lineshape
        which is a function of energy list and energy transfer.'''
        # make angular momentum ket object
        #eigenkets = [Ket(ei) for ei in self.eigenvectors]

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
        '''1D neutron spectrum without the Kf/Ki correction'''
        # make angular momentum ket object
        # eigenkets = [Ket(ei) for ei in self.eigenvectors]
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
        '''1D neutron spectrum without the Kf/Ki correction.
        LineshapeFunc must be a function with arguments of energy list and 
        central energy.'''
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
        intensity1D = self.neutronSpectrum(Earray, Temp, Ei, ResFunc,  gamma)

        # Scale by Debye-Waller Factor
        DWF = np.exp(1./3. * Qarray**2 * DebyeWaller**2)
        # Scale by form factor
        FormFactor = RE_FormFactor(Qarray,Ion)
        return np.outer(intensity1D, DWF*FormFactor)

    def normalizedNeutronSpectrum2D(self, Earray, Qarray, Temp, ResFunc, gamma, Ion, DebyeWaller=0):
        intensity1D = self.normalizedNeutronSpectrum(Earray, Temp, ResFunc,  gamma)

        # Scale by Debye-Waller Factor
        DWF = np.exp(1./3. * Qarray**2 * DebyeWaller**2)
        # Scale by form factor
        FormFactor = RE_FormFactor(Qarray,Ion)
        return np.outer(intensity1D, DWF*FormFactor)


    def _transition(self,ket1,ket2):  ## Correct, but slow.
        """Computes \sum_a |<|J_a|>|^2"""
        # Jx = Operator.Jx(ket1.j)
        # Jy = Operator.Jy(ket1.j)
        # Jz = Operator.Jz(ket1.j)
        # ax = np.dot(ket1.ket,np.dot(Jx.O,ket2.ket)) * np.dot(ket2.ket,np.dot(Jx.O,ket1.ket))
        # ay = np.dot(ket1.ket,np.dot(Jy.O,ket2.ket)) * np.dot(ket2.ket,np.dot(Jy.O,ket1.ket))
        # az = np.dot(ket1.ket,np.dot(Jz.O,ket2.ket)) * np.dot(ket2.ket,np.dot(Jz.O,ket1.ket))

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
        return 1/np.pi * (0.5*gamma)/((x-x0)**2 + (0.5*gamma)**2)

    def _voigt(self, x, x0, alpha, gamma):
        """ Return the Voigt line shape at x with Lorentzian component FWHM gamma
        and Gaussian component FWHM alpha."""
        sigma = (0.5*alpha) / np.sqrt(2 * np.log(2))
        return np.real(wofz(((x-x0) + 1j*(0.5*gamma))/sigma/np.sqrt(2))) / sigma\
                                                            /np.sqrt(2*np.pi)

    def _Re(self,value):
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
        '''prints eigenvectors and eigenvalues in a matrix'''
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
        '''prints eigenvectors and eigenvalues in the output that Latex can read'''
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
        """Prints <J_x>, <J_y>, and <J_z> for the ground state"""
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
        '''field should be a 3-component vector. Temps may be an array.
        Returns a three-component vector [M_x, M_y, M_z].
        Field should be in units of Tesla, and magnetization is calculated in Bohr Magnetons'''
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
        '''Computes susceptibility numerically with a numerical derivative.
        deltaField needs to be a scalar value.
        Returns a powder average value if Field is a scalar, and returns
        [Chi_x, Chi_y, Chi_z] if Field is a vector.
        Field should be in Tesla, and susceptibility is calculated in Bohr Magnetons
        per Tesla.'''
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
        # Compute susceptibility from perturbation theory, using MesotFurer eq. 11
        gJ = LandeGFactor(ion)
        muB = 5.7883818012e-2  # meV/T
        k_B = 8.6173303e-2  # meV/K

        # In this case, we assume powder average.

        # Jx = Operator.Jx(self.J)
        # Jy = Operator.Jy(self.J)
        # Jz = Operator.Jz(self.J)


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
        '''Returns g tensor computed numerically'''

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
        # Jx = Operator.Jx(self.J).O
        # Jy = Operator.Jy(self.J).O
        # Jz = Operator.Jz(self.J).O
        Jx = self.opttran.Jx
        Jy = self.opttran.Jy*1j
        Jz = self.opttran.Jz
        #print(vv1,'\n',vv2)
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
         '''Returns g tensor computed numerically from zeeman splitting'''
         Jx = Operator.Jx(self.J)
         Jy = Operator.Jy(self.J)
         Jz = Operator.Jz(self.J)

         #print(Jx)
         #print(Jy)
         muB = 5.7883818012e-2  # meV/T
         #mu0 = np.pi*4e-7       # T*m/A

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
        '''fits data to CEF parameters'''

        # define parameters
        # if len(self.B) != len(kwargs['coeff']):
        #     raise ValueError('coeff needs to have the same length as self.B')

        # Define function to be fit
        fun, p0, resfunc = makeFitFunction(chisqfunc, fitargs, **dict(kwargs, CFLevelsObject=self) )

        ############## Fit, using error function  #####################
        p_best = optimize.minimize(fun, p0, method=method)
        ###############################################################

        #print(fun(p_best.x))
        #print(chisqfunc(self, **kwargs))
        initialChisq, finalChisq = chisqfunc(self, **kwargs), fun(p_best.x)
        print('\rInitial err =', initialChisq, '\tFinal err =', finalChisq)
        
        result = resfunc(p_best.x)
        #print '\nFinal values: ', result
        result['Chisq'] = finalChisq
        return result

    def fitdata_GlobalOpt(self, chisqfunc, fitargs, **kwargs):
        '''fits data to CEF parameters using the basin hopping algorithm'''

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
        #print '\nFinal values: ', result
        result['Chisq'] = finalChisq
        return result


    def testEigenvectors(self):
        """Tests if eigenvectors are really eigenvectors"""
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

