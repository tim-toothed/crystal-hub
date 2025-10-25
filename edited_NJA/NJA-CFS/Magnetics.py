import scipy
import numpy as np


#======================= MAGNETIC PROPERTIES ==============================

class Magnetics():
    """
    Class for magnetic properties calculations.

    References:
    1. R. Boca, "theoretical foundations of molecular magnetism" 1999.

    Attributes:
    result (np.array): eigenvalue and eigenvectors of the Hamiltonian matrix.
    par_dict (dict): dictionary of parameters used in the calculation.
    calc (Calculation): instance of the Calculation class.
    basis (np.array): basis of the Hamiltonian matrix.
    Hterms (list): The contributions to the Hamiltonian matrix among: 
                    'Hee'=interelectronic repulsion, 'Hso'=spin-orbit coupling, 
                    'Hcf'=crystal-field/ligand-field contribution, 'Hz'=Zeeman splitting.
    kB (float): Boltzmann constant.
    mu0 (float): vacuum permeability.
    muB (float): Bohr magneton.
    """

    def __init__(self, calc, contributes, par, wordy=False):
        """
        Initializes the Magnetics class, computing the 0 field results.

        Parameters:
        calc (Calculation): instance of the Calculation class.
        contributes (list): The contributions to the Hamiltonian matrix among: 
                            'Hee'=interelectronic repulsion, 'Hso'=spin-orbit coupling,
                            'Hcf'=crystal-field/ligand-field contribution, 'Hz'=Zeeman splitting.
        par (dict): dictionary of parameters used in the calculation.
        """

        self.result = calc.MatrixH(contributes, **par, wordy=wordy)
        self.par_dict = calc.par_dict

        self.calc = calc
        self.basis = self.calc.basis
        self.Hterms = contributes

        self.kB = scipy.constants.k
        self.mu0 = scipy.constants.mu_0
        self.muB = scipy.constants.physical_constants['Bohr magneton'][0]/1.9865e-23

    def mag_moment(self, k=1, evaluation=True):
        #costruction of magnetic moment matrix as kL+geS
        #y component is divided by i (imaginary unit)
        """
        Computes the magnetic moment matrix.

        Parameters:
        k (int, optional): orbital reduction factor.
        evaluation (bool, optional): A flag indicating whether the parameters are symbolic or numerical. Default is True.
        
        Returns:
        matrix (numpy.ndarray): the magnetic moment matrix (complex).
        """

        if not evaluation:
            matrix = np.zeros((3, self.basis.shape[0],self.basis.shape[0]),dtype=object)
        else:
            matrix = np.zeros((3, self.basis.shape[0],self.basis.shape[0]),dtype='complex128')
        for i in range(self.basis.shape[0]):
            statei = self.basis[i]
            Si = statei[0]/2.
            Li = statei[1]
            Ji = statei[2]/2.
            MJi = statei[3]/2.
            seni = statei[4]
            labeli = self.calc.dic_LS[':'.join([f'{qq}' for qq in statei])]
            for j in range(0,i+1):
                statej = self.basis[j]
                Sj = statej[0]/2.
                Lj = statej[1]
                Jj = statej[2]/2.
                MJj = statej[3]/2.
                senj = statej[4]
                labelj = self.calc.dic_LS[':'.join([f'{qq}' for qq in statej])]
                H = Hamiltonian([seni,Li,Si,senj,Lj,Sj,Ji,MJi,Jj,MJj], [labeli,labelj], self.calc.conf, self.calc.dic_cfp, self.calc.tables, self.calc.dic_LS, self.calc.dic_LS_almost)  #self.conf è quella del main
                if Li==Lj and Si==Sj and seni==senj:

                    kL, gS = H.Zeeman(k=k, evaluation=evaluation, MM=True)

                    # from -1,0,+1 to x,y,z
                    matrix[0,i,j] += (kL[0]+gS[0] - (kL[2]+gS[2]))*1/(np.sqrt(2))
                    matrix[1,i,j] += (kL[0]+gS[0] + kL[2]+gS[2])*1j/(np.sqrt(2))
                    matrix[2,i,j] += kL[1]+gS[1]

                    if i!=j:
                        for kk in range(3):
                            matrix[kk,j,i] += np.conj(matrix[kk,i,j])

        matrix = - matrix #the minus sign is because mu = - kL - 2S

        if evaluation:
            matrix = np.round(matrix, 16)

        return matrix

    def calc_LS(self, k=1):
        #costruction of Lx,y,z and Sx,y,z operators
        """
        Computes the magnetic moment matrix.

        Parameters:
        k (int, optional): orbital reduction factor.

        Returns:
        L (numpy.ndarray): the orbital angular momentum matrix [x,y,z] (multiplied by k).
        S (numpy.ndarray): the spin angular momentum matrix [x,y,z] (multiplied by ge).
        """

        L = np.zeros((3, self.basis.shape[0],self.basis.shape[0]),dtype='complex128')
        S = np.zeros((3, self.basis.shape[0],self.basis.shape[0]),dtype='complex128')
        for i in range(self.basis.shape[0]):
            statei = self.basis[i]
            Si = statei[0]/2.
            Li = statei[1]
            Ji = statei[2]/2.
            MJi = statei[3]/2.
            seni = statei[4]
            labeli = self.calc.dic_LS[':'.join([f'{qq}' for qq in statei])]
            for j in range(0,i+1):
                statej = self.basis[j]
                Sj = statej[0]/2.
                Lj = statej[1]
                Jj = statej[2]/2.
                MJj = statej[3]/2.
                senj = statej[4]
                labelj = self.calc.dic_LS[':'.join([f'{qq}' for qq in statej])]
                H = Hamiltonian([seni,Li,Si,senj,Lj,Sj,Ji,MJi,Jj,MJj], [labeli,labelj], self.calc.conf, self.calc.dic_cfp, self.calc.tables, self.calc.dic_LS, self.calc.dic_LS_almost)  #self.conf è quella del main
                if Li==Lj and Si==Sj and seni==senj:

                    kL, gS = H.Zeeman(k=k, evaluation=True, MM=True)

                    # from -1,0,+1 to x,y,z
                    L[0,i,j] += (kL[0] - kL[2])*1/(np.sqrt(2))
                    L[1,i,j] += (kL[0] + kL[2])*1j/(np.sqrt(2))
                    L[2,i,j] += kL[1]

                    S[0,i,j] += (gS[0] - gS[2])*1/(np.sqrt(2))
                    S[1,i,j] += (gS[0] + gS[2])*1j/(np.sqrt(2))
                    S[2,i,j] += gS[1]

                    if i!=j:
                        for kk in range(3):
                            L[kk,j,i] += np.conj(L[kk,i,j])
                            S[kk,j,i] += np.conj(S[kk,i,j])

        return L, S

    #@staticmethod
    def effGval(self, levels, v_matrix=None):
        """
        Computes the effective g-matrix for the Kramers doublets specified.

        Parameters:
        levels (list): list of the Kramers doublet levels, e.g. [(1,2), (3,4)].
        v_matrix (numpy.ndarray, optional): eigenvectors of the Hamiltonian matrix. Default is None.

        Returns:
        np.sqrt(2*w) (numpy.ndarray): eigenvalues of the g-matrix.
        v (numpy.ndarray): eigenvectors of the g-matrix.
        """ 

        if v_matrix is None:
            par = self.par_dict
            matrix = self.calc.build_matrix(self.Hterms, **par) 
            result = from_matrix_to_result(matrix)
            v_matrix = np.copy(result[1:,:])
        else:
            pass

        levels = np.array(levels)

        mu_matrix = self.mag_moment(1)  
                                                             
        gexs = np.zeros(10, dtype='int32')
        ngexs = int((levels[1]-levels[0]+1)/2)

        for i in range(ngexs):
            gexs[2*i]=int(levels[0])+2*i-1
            gexs[2*i+1]=int(levels[0])+2*i

        G2 = np.zeros((3,3), dtype='complex128')
        for i in range(ngexs):
            j = gexs[2*i]
            idx1 = [j,j,j+1,j+1]
            idx2 = [j,j+1,j,j+1]
            gk = np.zeros((3,4,4), dtype='complex128')
            for ii in idx1:
                for jj in idx2:
                    for kk in range(3):
                        gk[kk,ii,jj] = np.dot(np.conj(v_matrix[:,ii]).T, np.dot(mu_matrix[kk,...],v_matrix[:,jj]))

            gx11 = gk[0,j,j]; gx12 = gk[0,j,j+1]; gx21 = gk[0,j+1,j]; gx22 = gk[0,j+1,j+1];
            gy11 = gk[1,j,j]; gy12 = gk[1,j,j+1]; gy21 = gk[1,j+1,j]; gy22 = gk[1,j+1,j+1];
            gz11 = gk[2,j,j]; gz12 = gk[2,j,j+1]; gz21 = gk[2,j+1,j]; gz22 = gk[2,j+1,j+1];

            G2[0,0]=gx11*gx11 + gx12*gx21 + gx21*gx12 + gx22*gx22
            G2[0,1]=gx11*gy11 + gx12*gy21 + gx21*gy12 + gx22*gy22
            G2[1,0]=G2[0,1]
            G2[0,2]=gx11*gz11 + gx12*gz21 + gx21*gz12 + gx22*gz22
            G2[2,0]=G2[0,2]
            G2[1,1]=gy11*gy11 + gy12*gy21 + gy21*gy12 + gy22*gy22
            G2[1,2]=gy11*gz11 + gy12*gz21 + gy21*gz12 + gy22*gz22
            G2[2,1]=G2[1,2]
            G2[2,2]=gz11*gz11 + gz12*gz21 + gz21*gz12 + gz22*gz22

        w,v = np.linalg.eigh(G2)

        return np.sqrt(2*w),v


    def susceptibility_field(self, fields, temp, delta=0.001, wordy=False):
        """
        Computes the scalar magnetization and magnetic susceptibility fields, for a set of field vectors (evenly sampled on a sphere) at a certain temperature.

        Parameters:
        fields (numpy.ndarray): field vectors (shape = N x 3, in T).
        temp (float): temperature in K.
        delta (float, optional): differentiation step in T. Default is 0.001.
        wordy (bool, optional): A flag indicating whether to print the results. Default is False.

        Returns:
        M_list (numpy.ndarray): scalar magnetization field (SI unit).
        susc_list (numpy.ndarray): scalar magnetic susceptibility field (SI unit).
        """

        par = self.par_dict

        M_list = np.zeros(fields.shape[0])
        susc_list = np.zeros(fields.shape[0])
        mu = np.zeros((fields.shape[0], self.basis.shape[0], 3), dtype='complex128')

        mu_matrix = self.mag_moment(1)

        for i in range(fields.shape[0]): 
            print(str(i)+' FIELD: '+str(fields[i])+'\n' if wordy else "", end = "")

            par['field'] = fields[i] 
            weights = fields[i]/np.linalg.norm(fields[i])
            matrix = self.calc.build_matrix(self.Hterms, **par) 
            result = from_matrix_to_result(matrix)
            E = (result[0,:].real-min(result[0,:].real)) 

            den = 0
            num = 0
            for ii in range(len(E)):
                for kk in range(3):
                    mu_single = np.dot(np.conj(result[1:,ii]).T, np.dot(mu_matrix[kk,...],result[1:,ii]))
                    if np.abs(np.copy(mu_single).imag)<1e-9:
                        mu[i,ii,0] += np.copy(mu_single.real)*weights[kk]
                    else:
                        print('complex')

                num += np.real(mu[i,ii,0]*np.exp(-E[ii]/(self.kB/1.9865e-23*temp)))
                den += np.real(np.exp(-E[ii]/(self.kB/1.9865e-23*temp)))
            E_av = num/den/3

            B_inc = fields[i]/np.linalg.norm(fields[i])*delta
            par['field'] = fields[i]+B_inc
            matrix = self.calc.build_matrix(self.Hterms, **par)  
            result = from_matrix_to_result(matrix)
            E = (result[0,:].real-min(result[0,:].real)) 
            den = 0
            num = 0
            for ii in range(len(E)):
                for kk in range(3):
                    mu_single = np.dot(np.conj(result[1:,ii]).T, np.dot(mu_matrix[kk,...],result[1:,ii]))
                    if np.abs(np.copy(mu_single).imag)<1e-9:
                        mu[i,ii,1] += np.copy(mu_single.real)*weights[kk]
                    else:
                        print('complex')

                num += np.real(mu[i,ii,1]*np.exp(-E[ii]/(self.kB/1.9865e-23*temp)))
                den += np.real(np.exp(-E[ii]/(self.kB/1.9865e-23*temp)))
            E_av_inc = num/den/3

            B_inc = fields[i]/np.linalg.norm(fields[i])*delta
            par['field'] = fields[i]-B_inc
            matrix = self.calc.build_matrix(self.Hterms, **par)  
            result = from_matrix_to_result(matrix)
            E = (result[0,:].real-min(result[0,:].real)) 
            den = 0
            num = 0
            for ii in range(len(E)):
                for kk in range(3):
                    mu_single = np.dot(np.conj(result[1:,ii]).T, np.dot(mu_matrix[kk,...],result[1:,ii]))
                    if np.abs(np.copy(mu_single).imag)<1e-9:
                        mu[i,ii,2] += np.copy(mu_single.real)*weights[kk]
                    else:
                        print('complex')

                num += np.real(mu[i,ii,2]*np.exp(-E[ii]/(self.kB/1.9865e-23*temp)))
                den += np.real(np.exp(-E[ii]/(self.kB/1.9865e-23*temp)))
            E_av_incm = num/den/3

            M_list[i] = E_av*self.muB*1.9865e-23
            susc_list[i] = (E_av_inc - E_av_incm)/(2*delta)*self.mu0*self.muB*1.9865e-23

            print('M '+str(M_list[i])+'\n' if wordy else "", end = "")
            print('chi '+str(susc_list[i])+'\n' if wordy else "", end = "")

        return M_list, susc_list
    

    def susceptibility_B_copy(self, fields, temp, delta=0.001, wordy=False):
        """
        Computes the magnetic susceptibility tensor as the derivative of the magnetization vector in x,y,z,
        for a set of field vectors (evenly sampled on a sphere). The values are then averaged among the 
        different directions of the magnetic field vector, hence it is generally enough to pass a single field vector 
        as e.g. np.array([[0,0,1]]) (fields needs to be bidimensional).

        Parameters:
        fields (numpy.ndarray): field vectors (shape = N x 3, in T). 
        temp (float): temperature in K.
        delta (float, optional): differentiation step in T. Default is 0.001.
        wordy (bool, optional): A flag indicating whether to print the results. Default is False.

        Returns:
        chi_tensor (numpy.ndarray): magnetic susceptibility tensor (SI unit).
        Mav (numpy.ndarray): magnetization vector (SI unit).
        """

        def M_vector_in(field_vec):
            mu = np.zeros((self.basis.shape[0], 3), dtype='complex128')
            par['field'] = field_vec
            matrix = self.calc.build_matrix(self.Hterms, **par)  #hamiltonian matrix for a field vector B
            result = from_matrix_to_result(matrix)
            E = (result[0,:].real-min(result[0,:].real)) #* 1.9865e-23
            E -= min(E)

            M = np.zeros(3)

            for kk in range(3):
                den = 0
                num = 0
                for ij in range(len(E)):
                    mu_single = np.dot(np.conj(result[1:,ij]).T, np.dot(mu_matrix[kk,...],result[1:,ij]))   # <i|mu_kk|i> for i in range N and for kk=x,y,z
                    if np.abs(np.copy(mu[ij,kk]).imag)<1e-15:
                        mu[ij,kk] += np.copy(mu_single.real)
                    else:
                        print('complex',mu_single)  #just to check that the values are real and everything is okay

                    num += np.real(mu[ij,kk]*np.exp(-E[ij]/(self.kB/1.9865e-23*temp)))
                    den += np.real(np.exp(-E[ij]/(self.kB/1.9865e-23*temp)))

                M[kk] = num/den

            M = np.round(M, 16)
            return M

        par = self.par_dict
        try:
            k=par['k']
        except:
            k=1

        Mag_vector = np.zeros((fields.shape[0],3))
        chi = np.zeros((fields.shape[0], 3, 3))
        mu_matrix = self.mag_moment(k)  #computes the magnetic moments matrix <i|kL + geS|j>

        for i in range(fields.shape[0]):
            print('FIELD: '+str(fields[i])+'\n' if wordy else "", end = "")

            M = M_vector_in(fields[i])

            Mag_vector[i,:] = M
            field_inc = np.zeros((3,3))    #matrix of fields increment, one incremented field vector for the three components: x,y,z
            field_dec = np.zeros_like(field_inc)
            M_inc = np.zeros_like(field_inc, dtype='float64')
            M_dec = np.zeros_like(field_inc, dtype='float64')
            for comp in range(3):
                field_inc[comp] = np.copy(fields[i])
                field_inc[comp, comp] += delta
                M_inc[comp,:] = M_vector_in(np.round(field_inc[comp],16))   #computation of the magnetization vector (see above) for every magnetic field increment
                field_dec[comp] = np.copy(fields[i])
                field_dec[comp, comp] -= delta
                M_dec[comp,:] = M_vector_in(np.round(field_dec[comp],16))   #same for decrement

            for kki in range(3):
                for kkj in range(3):
                    chi[i,kki,kkj] = ((M_inc[kkj,kki]-M_dec[kkj,kki])/(2*delta))*self.mu0*self.muB*1.9865e-23

            print('M (BM)\n'+str(M)+'\n' if wordy else "", end = "")
            print('M_inc (BM)\n'+str(M_inc)+'\n' if wordy else "", end = "")
            print('chi (m3)\n'+str(chi[i,:,:])+'\n' if wordy else "", end = "")

        chi_tensor = np.zeros((3,3))
        chi_tensor = np.sum(chi, axis=0)/fields.shape[0]   #since I've calculated the tensor for different directions of the magnetic field vector, here there is the average calculation
        Mav = np.sum(Mag_vector, axis=0)/fields.shape[0]    #same for magnetization

        return chi_tensor, Mav*self.muB*1.9865e-23

    # @staticmethod
    # @cron
    ### COPY OF THE NUMBA FUNCTIONS ###
    def dfridr(self, func, x, h, idxi, shape, fargs):
        """
        returns the derivative of the function at a point x by Ridders' method of polynomial extrapolation. 
        The value h is input as an estimated initial stepsize. It has to be bigger than the one you would pass to a standard numerical differentiation method. 
        The stepsize is decreased by CON at each iteeration. Max size of tableau is set by NTAB.
        See NJA-CFS documentation for references.

        Parameters:
        func (function): function to differentiate.
        x (numpy.ndarray): differentiation variable.
        h (float): step size.
        idxi (int): index of the differentiation variable.
        shape (tuple): shape of the output.
        fargs (tuple): additional arguments for the function.

        Returns:
        result (numpy.ndarray): result of the differentiation.
        err (float): error estimation of the differentiation.
        """

        CON = h*2 #* 2  #10  #if this is too high the error at the end will be higher, but if it's too low the result will be always 0
        CON2 = CON * CON
        NTAB = 10  #10
        SAFE = 2  #2
        a = np.zeros((NTAB, NTAB)+shape[1:])

        hh = h
        zero = 1e-16

        dx = np.copy(x)
        dx[idxi] += hh
        sx = np.copy(x)
        sx[idxi] -= hh
        if 2*hh!=0:
            a[0,0,...] = ((func(dx,*fargs)-func(sx,*fargs))/(2*hh))
        else:
            a[0,0,...] = ((func(dx,*fargs)-func(sx,*fargs))/zero)

        err = np.inf
        result = None

        for i in range(1, NTAB):
            hh /= CON
            dx = np.copy(x)
            dx[idxi] += hh
            sx = np.copy(x)
            sx[idxi] -= hh
            if 2*hh!=0:
                a[0,i,...] = ((func(dx,*fargs)-func(sx,*fargs))/(2*hh))
            else:
                a[0,i,...] = ((func(dx,*fargs)-func(sx,*fargs))/zero)
            fac = CON2
            for j in range(1, i):
                if (fac - 1)!=0:
                    a[j, i,...] = (a[j - 1, i,...] * fac - a[j - 1, i - 1,...]) / (fac - 1)
                else:
                    a[j, i,...] = (a[j - 1, i,...] * fac - a[j - 1, i - 1,...]) / zero
                fac *= CON2
                errt = max(norm(a[j, i,...] - a[j - 1, i,...]), norm(a[j, i,...] - a[j - 1, i - 1,...]))
                if errt <= err:
                    err = errt
                    result = a[j, i,...]
            if norm(a[i, i,...] - a[i - 1, i - 1,...]) >= SAFE * err:
                return result, err

        return result, err

    #@cron
    ### COPY OF THE NUMBA FUNCTIONS ###
    def susceptibility_B_ord1(self, fields, temp, basis, LF_matrix, delta=1.):
        """
        Computes the magnetic susceptibility tensor as the derivative of the magnetization vector in x,y,z,
        for a set of field vectors (evenly sampled on a sphere). The values are then averaged among the 
        different directions of the magnetic field vector, hence it is generally enough to pass a single field vector
        as e.g. np.array([[0,0,1]]) (fields needs to be bidimensional).
        The differentiation method is based on Ridders' method of polynomial extrapolation. 
        The value h is input as an estimated initial stepsize (it has to be larger that the one you would pass to a standard differentiation procedure).
        An estimate of the error is also computed.

        Parameters:
        fields (numpy.ndarray): field vectors (shape = N x 3, in T).
        temp (float): temperature in K.
        basis (numpy.ndarray): basis of the Hamiltonian matrix.
        LF_matrix (numpy.ndarray): 0 field Hamiltonian matrix.
        delta (float, optional): differentiation step in T. Default is 1.

        Returns:
        chi_tensor (numpy.ndarray): magnetic susceptibility tensor (SI unit).
        err_tensor (numpy.ndarray): error estimation for the computation of the magnetic susceptibility tensor (SI unit).
        """

        mu0 = 1.25663706212e-06
        muB = 0.4668517532494337
        Jconv = 1.9865e-23

        mu_matrix = mag_moment(basis)  #complex128[:,:,:]
        chi = np.zeros((fields.shape[0], 3, 3), dtype='float64')
        err = np.zeros_like(chi)
        if len(fields.shape)<2:
            fields = np.array([fields])
        for i in range(fields.shape[0]):
            for idx in range(3):
                chi_comp, err_comp = self.dfridr(M_vector, fields[i], delta, idx, chi.shape[1:], fargs=(mu_matrix, LF_matrix, basis, temp))
                chi[i,idx] = chi_comp * mu0*muB*Jconv
                err[i,idx] = np.ones(chi_comp.shape)*err_comp * mu0*muB*Jconv

        chi_tensor = np.zeros((3,3))
        chi_tensor = np.sum(chi, axis=0)/fields.shape[0]
        err_tensor = np.sum(err, axis=0)/fields.shape[0]

        return chi_tensor, err_tensor
