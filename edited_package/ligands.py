import numpy as np
from scipy import optimize

from edited_package.lattice_class import lattice
from edited_package.plot_ligands import exportLigandCif
from edited_package.constants import calculate_tesseral_harmonic, theta, calculate_radial_integral_RE, Constant, LStheta, PFalpha, PFbeta, calculate_radial_integral_TM, ION_NUMS_RARE_EARTH, is_half_filled
from edited_package.stevens_operators import StevensOp, LS_StevensOp
from edited_package.create_fit_function import makeFitFunction
from edited_package.cf_levels import CFLevels, LS_CFLevels
from edited_package.operators import LSOperator


class Ligands:
    """For doing point-charge calculations"""
    def __init__(self,ion,ligandPos, latticeParams=None, ionPos=[0,0,0]):
        """Creates array of ligand bonds in cartesian coordinates"""
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
        '''rotates the ligand bonds so that the new axis is in the direction of the old axis'''
        rotationAxis = np.cross(newaxis,oldaxis)
        rotationAngle = np.arccos(np.dot(newaxis,oldaxis)/(np.linalg.norm(newaxis)*np.linalg.norm(oldaxis)))
        self.bonds = np.array([self._rotateMatrix(b,rotationAxis,rotationAngle) for b in self.bonds])

    def rotateLigandsZ(self, oldaxis):
        '''rotates the ligand bonds around the z axis so that oldaxis 
        becomes the x axis'''
        zrotation = np.arctan(oldaxis[1]/oldaxis[0])
        self.bonds = np.array([self._rotateMatrix(b,np.array([0,0,1]),-zrotation) for b in self.bonds])


    def _rotateMatrix(self,matrixin,axis,angle):
        """taken from http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/"""
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


    def exportCif(self, filename):
        exportLigandCif(self, filename)


    def PointChargeModel(self, symequiv=None, LigandCharge=-2,IonCharge=3, printB = True, 
                            suppressminusm = False, ionL=None):
        '''Create point charge model of the crystal fields of a rare-earth ion.
        Returns a CFLevels object with the hamiltonian defined.
        Define LigandCharge in units of e.'''

        self.IonCharge = IonCharge
        # Lock suppressmm into whatever it was when PointChargeModel was first called.
        try: self.suppressmm
        except AttributeError:
            self.suppressmm = suppressminusm


        if symequiv == None:
            # charge = IonCharge*[LigandCharge]*len(self.bonds)
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
                #charge[i] = IonCharge*LigandCharge[se]
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
        #for n,m in [[n,m] for n in range(2,8,2) for m in range(-n,n+1)]:
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


    def FitChargesNeutrons(self, chisqfunc, fitargs, method='Powell', **kwargs):
        '''This is the old name. I keep it around so that the original code 
        will run with it.'''
        return self.FitCharges(chisqfunc, fitargs, method='Powell', **kwargs)

    def FitCharges(self, chisqfunc, fitargs, method='Powell', **kwargs):
        '''fits data'''

        # Define function to be fit
        fun, p0, resfunc = makeFitFunction(chisqfunc, fitargs, **dict(kwargs, LigandsObject=self) )

        print('\tFitting...')
        ############## Fit, using error function  #####################
        p_best = optimize.minimize(fun, p0, method=method)
        ###############################################################

        try:
            initialChisq, finalChisq = fun(p0), fun(p_best.x)
            finalvals = resfunc(p_best.x)
        except IndexError:
            initialChisq, finalChisq = fun(p0), fun([float(p_best.x)])
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

###################################
### Same class, but in the LS basis
###################################

# LS basis functions are a set of quantum mechanical wave functions
# that describe the electronic states of an atom 
# in the L-S (Russell-Saunders) coupling scheme. 
# This model is used for lighter atoms where 
# the electrostatic repulsion between electrons
# is stronger than the spin-orbit interaction.

class LS_Ligands:
    """For doing point-charge calculations in LS basis"""
    def __init__(self,ion, ligandPos, SpinOrbitCoupling, latticeParams=None, ionPos=[0,0,0]):
        """Creates array of ligand bonds in cartesian coordinates.
        'ion' can either be the name of the ion or a list specifying L and S.
        For example, it could be 'Ni3+', or ['Ni3+', 0.5, 1]"""
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
        '''rotates the ligand bonds so that the new axis is in the direction of the old axis'''
        rotationAxis = np.cross(newaxis,oldaxis)
        rotationAngle = np.arccos(np.dot(newaxis,oldaxis)/(np.linalg.norm(newaxis)*np.linalg.norm(oldaxis)))
        self.bonds = np.array([self._rotateMatrix(b,rotationAxis,rotationAngle) for b in self.bonds])

    def rotateLigandsZ(self, oldaxis):
        '''rotates the ligand bonds around the z axis so that oldaxis 
        becomes the x axis'''
        zrotation = np.arctan(oldaxis[1]/oldaxis[0])
        self.bonds = np.array([self._rotateMatrix(b,np.array([0,0,1]),-zrotation) for b in self.bonds])


    def _rotateMatrix(self,matrixin,axis,angle):
        """taken from http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/"""
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
        '''Create point charge model of the crystal fields of a rare-earth ion.
        Returns a CFLevels object with the hamiltonian defined.
        Define LigandCharge in units of e.'''

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
        ''' For transition metals:
        Create point charge model of the crystal fields.
        Returns a CFLevels object with the hamiltonian defined.
        Define LigandCharge in units of e.'''
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
        ''' For transition metals if the radial integrals are not in PyCrystalField (d5 ions)
        Create point charge model of the crystal fields.
        Returns a CFLevels object with the hamiltonian defined.
        Define LigandCharge in units of e.'''

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
        '''fits neutron data'''

        # Define function to be fit
        fun, p0, resfunc = makeFitFunction(chisqfunc, fitargs, **dict(kwargs, LigandsObject=self) )

        print('\tFitting...')
        ############## Fit, using error function  #####################
        p_best = optimize.minimize(fun, p0, method=method)
        ###############################################################

        initialChisq, finalChisq = fun(p0), fun(p_best.x)

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
