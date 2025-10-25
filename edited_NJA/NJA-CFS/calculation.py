import scipy
import numpy as np
from pprint import pprint
import warnings
import crystdat

class calculation(): 
    """
    This is the main class of the code. Here the Hamiltonian matrix is computed and diagonalized.
    Additionally, in this class the basis set can be reduced and the wavefunctions optimized.

    Attributes:
    conf (str): The electron configuration.
    l (int): The one-electron orbital quantum number.
    n (int): The number of electrons in the almost-closed-shell configuration.
    N (int): The number of electrons in the configuration.
    closed (bool): A flag indicating whether the configuration is almost-closed-shell or not.
    basis (numpy.ndarray): The complete basis set for the configuration.
    dic_LS (dict): A dictionary for the LS-coupling scheme.
    basis_l (numpy.ndarray): The labels for the states.
    basis_l_JM (numpy.ndarray): The labels for the states (with MJ).   
    microst (int): The number of microstates.
    ground (bool): A flag indicating whether the basis set includes only the ground multiplet.
    dic_cfp (dict): A dictionary containing the coefficients of fractional parentage (cfp) for the configuration.
    tables (dict): A dictionary containing the calculated RMEs.
    dic_LS_almost (dict): An inverse dictionary for the LS-coupling scheme for the almost-closed-shell configuration.
    dic_ee (dict): A dictionary containing the electron-electron integrals. 
    """


    def __init__(self, conf, ground_only=False, TAB=False, wordy=True):
        """
        Initializes the calculation object through the definition of the type of configuration 
        (almost closed shell or not), the basis set and additional dictionaries.
        Additionally, it defines whether the RMEs are computed explicitly or read from tables 
        (the computation of the RMEs will result in an increase of the computational time).

        Parameters:
        conf (str): The electron configuration.
        ground_only (bool, optional): A flag indicating whether the basis set includes only the ground multiplet. Default is False.
        TAB (bool, optional): A flag indicating whether the RMEs are calculated or read from tables. Default is False.
        wordy (bool, optional): A flag indicating whether to print the output. Default is True.
        """

        conf_list_d = ['d1','d2','d3','d4','d5','d6','d7','d8','d9']
        conf_list_f = ['f1','f2','f3','f4','f5','f6','f7','f8','f9','f10','f11','f12','f13']
        if conf not in conf_list_d and conf not in conf_list_f:
            raise ValueError('Configuration not valid')

        if conf[0]=='d':
            self.l = 2
        else:
            self.l = 3

        self.conf = conf  #questa deve rimanere sempre quella del main
        self.n = int(conf[1:])  #elettroni closed conf
        self.N = int(conf[1:])  #numero di elettroni vero
        self.closed = False
        if self.l == 2 and self.n>5:
            self.closed = True
            self.n = almost_closed_shells(conf)
            stringa = 'd'+str(self.n+1)
        elif self.l == 2 and self.n<=5:
            stringa = conf
        elif self.l == 3 and self.n>7:
            self.closed = True
            self.n = almost_closed_shells(conf)
            stringa = 'f'+str(self.n+1)
        elif self.l == 3 and self.n<=7:
            stringa = conf

        self.basis, self.dic_LS, self.basis_l, self.basis_l_JM = Full_basis(conf)

        if ground_only==True:
            self.ground = True
            self.basis, self.dic_LS, self.basis_l, self.basis_l_JM = self.ground_state_calc()
        else:
            self.ground = False

        print('\nConfiguration: '+conf+'\n' if wordy else "", end = "")
        self.microst = int(scipy.special.binom(2*(self.l*2+1), self.N))
        print('Number of microstates: '+str(self.microst)+'\n' if wordy else "", end = "")
        if self.closed == True:
            print('Almost closed Shell, corresponding configuration: '+conf[0]+str(self.n)+'\n' if wordy else "", end = "")
        if ground_only==True:
            print('Ground state only calculation\nBasis set reduced to: '+str(self.basis.shape[0])+'\n' if wordy else "", end = "")

        self.dic_cfp, self.tables, self.dic_LS_almost, self.dic_ee = self.Tables(TAB, stringa, conf)

    def ground_state_calc(self, ground=None):
        #basis --> array: n. microstates x [2S, L, 2J, 2M, sen (,count)]
        #dic_LS --> dict: '2S:L:2J:2M:sen:count': label as N. and K.
        #basis_l --> list: 2S+1 L (J)
        #basis_l_JM --> list: 2S+1 L (J) MJ
        """
        Reduces the basis set to the ground multiplet.

        Parameters:
        ground (str, optional): The ground multiplet. Default is None.

        Returns:
        basis_red (numpy.ndarray): The reduced basis set for the ground multiplet.
        dic_LS_red (dict): A dictionary for the LS-coupling scheme for the ground multiplet.
        basis_l_red (numpy.ndarray): The labels for the states for the ground multiplet (hence, only one type).
        basis_l_JM_red (numpy.ndarray): The labels for the states (with MJ) for the ground multiplet.
        """
        if ground is None:
            ground = ground_term_legend(self.conf) #conf NOT almost closed
        term_num = [(int(ground[0])-1), state_legend(ground[1]), eval(ground[ground.index('(')+1:ground.index(')')])*2]
        basis_red = []
        dic_LS_red = {}
        for i in range(self.basis.shape[0]):
            if list(self.basis[i,:3])==term_num:
                basis_red.append(list(self.basis[i]))
                lista = list(self.basis[i])
                dic_LS_red[':'.join(f'{n}' for n in lista)]=self.dic_LS[':'.join(f'{n}' for n in lista)]
        basis_red = np.array(basis_red)
        basis_l_red = [item for item in self.basis_l if item==ground]
        basis_l_JM_red = [item for item in self.basis_l_JM if item[:(item.index(')')+1)]==ground]

        return basis_red, dic_LS_red, basis_l_red, basis_l_JM_red

    def Tables(self, TAB, stringa, conf):
        """
        Reads the RMEs from tables or computes them explicitly.

        Parameters:
        TAB (bool): A flag indicating whether the RMEs are calculated or read from tables.
        stringa (str): The corresponding configuration.
        conf (str): The electron configuration.

        Returns:
        dic_cfp (dict): A dictionary containing the coefficients of fractional parentage (cfp) for the configuration.
        tables (dict): A dictionary containing the calculated RMEs.
        dic_LS_almost (dict): An inverse dictionary for the LS-coupling scheme for the almost-closed-shell configuration.
        dic_ee (dict): A dictionary containing the electron-electron integrals.
        """
        
        if conf=='d1' or conf=='f1':
            if conf=='d1':
                tables = {'U2': {'2D': {'2D': np.float64(1.0)}},
                        'U4': {'2D': {'2D': np.float64(1.0)}},
                        'V11': {'2D': {'2D': np.float64(1.224744871391589)}}}
                dic_ee = {'2D': {'2D': [0,0,0,0]}}
                dic_LS_almost = None
                dic_ee = None
                dic_cfp = None
            else:
                tables = {'U1': {'2F': {'2F': np.float64(1.0)}},
                    'U2': {'2F': {'2F': np.float64(1.0)}},
                    'U3': {'2F': {'2F': np.float64(1.0)}},
                    'U4': {'2F': {'2F': np.float64(1.0)}},
                    'U5': {'2F': {'2F': np.float64(1.0)}},
                    'U6': {'2F': {'2F': np.float64(1.0)}},
                    'V11': {'2F': {'2F': np.float64(1.224744871391589)}}}
                dic_ee = {'2F': {'2F': [0,0,0,0]}}
                dic_LS_almost = None
                dic_ee = None
                dic_cfp = None
        else:
            if TAB==False:
                dic_cfp = cfp_from_file(stringa)
                tables = None
                if self.closed==True:
                    dic_LS_almost = Full_basis(conf[0]+str(self.n+1))[1]   #, self.n+1)[1]  #questa è la dic_LS della configurazione corrispondente per i cfp
                else:
                    dic_LS_almost = None
            else:
                dic_cfp = None
                tables = read_matrix_from_file(self.conf, self.closed)
                dic_LS_almost = None
            dic_ee = None
            if self.l==3 and self.n>2:
                dic_ee = read_ee_int(self.conf, self.closed)

        return dic_cfp, tables, dic_LS_almost, dic_ee

    def reduce_basis(self, conf, roots, dic=None, contributes=None, wordy=False):
        # reduce the basis set according to spin multiplicity
        # for dN configurations the Hund's rules are applied
        # for fN configurations the weak field - High Spin situation is considered and 
        # default parameters are taken from Ma, C. G., Brik, M. G., Li, Q. X., & Tian, Y. (2014) Journal of alloys and compounds, 599, 93-101.

        # basis --> array: n. microstates x [2S, L, 2J, 2M, sen (,count)]
        # dic_LS --> dict: '2S:L:2J:2M:sen:count': label as N. and K.
        # basis_l --> list: 2S+1 L (J)
        # basis_l_JM --> list: 2S+1 L (J) MJ
        """
        Reduces the basis set according to the spin multiplicity.
        1. for dN configurations the Hund's rules are applied
        2. for fN configurations the weak field - High Spin situation is considered and 
        The default parameters are taken from Ma, C. G., Brik, M. G., Li, Q. X., & Tian, Y. (2014) Journal of alloys and compounds, 599, 93-101.

        Parameters:
        conf (str): The electron configuration.
        roots (list, optional): A list containing the selected multiplicities as list of (sum(2*L+1), 2*S+1) pairs.
        dic (dict, optional): A dictionary containing the parameters for the free ion. Default is None.
        contributes (list, optional): A list containing the contributions to the matrix. Default is None.
        """

        print('Performing basis reduction... \n' if wordy else "", end = "")

        def state_select(ground, basis, basis_l, basis_l_JM):

            term_num = [(int(ground[0])-1), state_legend(ground[1])]#, eval(ground[ground.index('(')+1:ground.index(')')])*2]
            basis_new = []
            basis_l_red = []
            basis_l_JM_red = []
            dic_LS_red = {}
            for ii in range(basis.shape[0]):
                if np.equal(basis[ii,:len(term_num)],term_num).all():
                    # print(basis[ii,:], term_num, basis_l[ii], basis_l_JM[ii])
                    if list(basis[ii]) not in basis_new:
                        basis_new.append(list(basis[ii]))
                        lista = list(basis[ii])
                        dic_LS_red[':'.join(f'{n}' for n in lista)]=self.dic_LS[':'.join(f'{n}' for n in lista)]
                        basis_l_red.append(basis_l[ii])
                        basis_l_JM_red.append(basis_l_JM[ii])
            basis_new = np.array(basis_new)
            
            basis_update = []
            basis_l_update = []
            basis_l_JM_update = []
            for ii in range(basis.shape[0]):
                if not any((basis[ii] == x).all() for x in basis_new):
                    basis_update.append(list(basis[ii]))
                    basis_l_update.append(basis_l[ii])
                    basis_l_JM_update.append(basis_l_JM[ii])
            basis_update = np.array(basis_update)

            return basis_update, basis_l_update, basis_l_JM_update, basis_new, dic_LS_red, basis_l_red, basis_l_JM_red

        # the basis must be complete
        if self.ground:
            raise ValueError("Basis set reduction not allowed in ground-only calculation on "+conf+" configuration")

        if conf[0]=='d':  # dN configurations follows Hund's rule
            basis_states = self.basis[:, :3]
            max_proj = np.unique(basis_states, axis=0)
            if int(conf[1:])>5:
                indices = np.lexsort((-max_proj[:, 2], -max_proj[:, 1], -max_proj[:, 0]))
            else:
                indices = np.lexsort((max_proj[:, 2], -max_proj[:, 1], -max_proj[:, 0]))  #the last one is the first criteria
            max_proj = max_proj[indices]
        else:   # fN configurations are in the weak field - High Spin situation

            if dic is None:
                dic = free_ion_param_f_HF(conf)
            if contributes is None:
                matrix = self.build_matrix(['Hee'], **dic)
            else:
                matrix = self.build_matrix(contributes, **dic)

            w,v = diagonalisation(matrix)
            v = np.round(v, 16)  #numbers are saved as complex128 data type
            result = np.vstack((w,v))
            result = np.copy(result[:, result[0,:].argsort()])

            projected = projection_basis(result[1:,:], self.basis_l)  #this gives just the order of spectroscopic terms of the free ion

            max_proj = []
            for i in projected.keys():
                s_proj = projected[i]
                dlabel = list(s_proj.keys())
                dvalue = np.array(list(s_proj.values()))
                max_proj.append(dlabel[dvalue.argmax()])

            unique_max_proj = []
            for item in max_proj:
                if item not in unique_max_proj:
                    unique_max_proj.append(item)
            max_proj = unique_max_proj

        basis_update = np.copy(self.basis)
        basis_l_update = np.copy(self.basis_l)
        basis_l_JM_update = np.copy(self.basis_l_JM)
        nroots = np.zeros(len(roots))
        for i in range(len(roots)):
            for j in range(len(max_proj)):
                if nroots[i] < roots[i][0]*roots[i][1] and int(max_proj[j][0])==roots[i][1]:
                    nroots[i] += eval(max_proj[j][0])*(2*state_legend(max_proj[j][1])+1) #(eval(max_proj[j][max_proj[j].index('(')+1:max_proj[j].index(')')]))*2 +1 #eval(max_proj[j][0])*(2*state_legend(max_proj[j][1])+1) #
                    basis_update, basis_l_update, basis_l_JM_update, basis_proj, dic_LS_proj, basis_l_proj, basis_l_JM_proj = state_select(max_proj[j], basis_update, basis_l_update, basis_l_JM_update)
                    
                    if basis_update.size > 0:
                        if i==0 and 'basis_red' not in locals():
                            basis_red = np.copy(basis_proj)
                            dic_LS_red = dic_LS_proj.copy()
                            basis_l_red = basis_l_proj.copy()
                            basis_l_JM_red = basis_l_JM_proj.copy()
                        else:
                            basis_red = np.vstack((basis_red, basis_proj))
                            dic_LS_red.update(dic_LS_proj)
                            basis_l_red += basis_l_proj
                            basis_l_JM_red += basis_l_JM_proj

        self.basis, self.dic_LS, self.basis_l, self.basis_l_JM = basis_red, dic_LS_red, basis_l_red, basis_l_JM_red

        print('Calculation on reduced set\nBasis set reduced to: '+str(self.basis.shape[0])+'\n' if wordy else "", end = "")


    #@njit
    #@cron
    def MatrixH(self, elem, F0=0, F2=0, F4=0, F6=0, zeta=0, k=1, dic_V=None, dic_bkq = None,dic_AOM = None, PCM = None, field = [0.,0.,0.], cfp_angles = None,wordy=False,
                      Orth=False, Norm=False, eig_opt=False, ground_proj=False, return_proj=False, save_label=False, save_LF=False, save_matrix=False):
        """
        This is the core of the calculation. This functions converts the different CFPs formalism in the Wybourne one, it optimizes the wavefunction
        composition, calls the function that builts the Hamiltonian and computes the wavefunctions composition.

        Parameters:
        elem (list): The contributions to the Hamiltonian matrix among: 
                    'Hee'=interelectronic repulsion, 'Hso'=spin-orbit coupling, 
                    'Hcf'=crystal-field/ligand-field contribution, 'Hz'=Zeeman splitting.
        F0 (float, optional): The F0 Slater-Condon parameter in cm^{-1}. Default is 0.
        F2 (float, optional): The F2 Slater-Condon parameter in cm^{-1}. Default is 0.
        F4 (float, optional): The F4 Slater-Condon parameter in cm^{-1}. Default is 0.
        F6 (float, optional): The F6 Slater-Condon parameter in cm^{-1}. Default is 0.
        zeta (float, optional): The spin-orbit coupling parameter in cm^{-1}. Default is 0.
        k (float, optional): The orbital reduction factor. Default is 1.
        dic_V (dict, optional): A dictionary containing one-electron ligand field matrix elements. Default is None. 
                                dic_V = {'11':0.0, '12':0, ...}, with key='ij', with 1<=i<=2*l+1 and j<=i.
        dic_bkq (dict, optional): A dictionary containing the crystal field parameters in Wybourne formalism. Default is None.
                                dic_bkq = {'k':{'q':0.0, '-q':0.0}, ...}, with k=2,4,..2*l and q=-k,-k+1,...k-1,k.
        dic_AOM (dict, optional): A dictionary containing the crystal field parameters in AOM formalism. Default is None.
                                dic_AOM = {'<ligand name>:[e_{sigma}, e_{pic}, e_{pis}, theta, phi, chi]}, where the angles are in degrees.
        PCM (list, optional): A list for the definition of the PCM for crystal field calculation. Default is None.
                              The first element of the list is the array dtype=object where each row is: '<ligand name>', x, y, z, charge.
                              The charge is expressed as fraction of electronic charge. The coordinates are in Angstrom. 
                              The metal center is assumed in (0,0,0). 
                              The second element in the list is a boolean that defines if the coordinates are expressed as cartesian (False) or spherical (True).
                              The third element is again a boolean and defines if the Sternheimer shielding parameters are used (True) or not (False).
        field (list, optional): The magnetic field vector in Tesla. Default is [0.,0.,0.].  
        cfp_angles (list, optional): A list containing the Euler angles or quaternions for the rotation of the crystal field parameters. Default is None.
        wordy (bool, optional): A flag indicating whether to print the output. Default is False.
        Orth (bool, optional): A flag indicating whether to perform the orthogonality check. Default is False.
        Norm (bool, optional): A flag indicating whether to perform the normalization check. Default is False.
        eig_opt (bool, optional): A flag indicating whether to optimize the wavefunction. Default is False.
        ground_proj (bool, optional): A flag indicating whether to wavefunctions composition is computed or not. Default is False.
        return_proj (bool, optional): A flag indicating whether to return the wavefunctions composition. Default is False.
                                      If the calculation is perfomed in the complete basis set or on a reduced one (different from just considering the ground state only with ground_only=True)
                                      the projection is a list of two dictionaries. The first without, the second with, the MJ label.
        save_label (bool, optional): A flag indicating whether to save the labels of the states. Default is False.
        save_LF (bool, optional): A flag indicating whether to save the Hamiltonian matrix (without the Zeeman contribution). Default is False.
        save_matrix (bool, optional): A flag indicating whether to save the Hamiltonian matrix. Default is False.

        Returns:
        result (numpy.ndarray): The eigenvalues and eigenvectors of the Hamiltonian matrix. The first row of the matrix are the eigenvalues.
                                Each column of the matrix is an eigenvector.
        projected (dict): The wavefunctions composition. The states are labeled as progressive numbers.
        """

        print('\nPerforming calculation with the following contributions: \n' if wordy else "", end = "")
        print(str(elem)+'\n' if wordy else "", end = "")
        if 'Hz' in elem:
            print('Magnetic field: \n'+f'{field[0]:.4e}'+' '+f'{field[1]:.4e}'+' '+f'{field[2]:.4e}'+' T\n' if wordy else "", end = "")

        #choice of conventions
        if 'Hcf' in elem:
            if dic_bkq is not None:
                pass
            elif dic_bkq is None and dic_V is not None:
                dic_bkq = from_Vint_to_Bkq(dic_V, self.conf)
            elif dic_bkq is None and dic_V is None and dic_AOM is not None:
                dic_V = from_AOM_to_Vint(dic_AOM, self.conf)
                dic_bkq = from_Vint_to_Bkq(dic_V, self.conf)
            elif dic_bkq is None and dic_V is None and dic_AOM is None and PCM is not None:
                dic_bkq = calc_Bkq(PCM[0], self.conf, PCM[1], PCM[2])
            else:
                print('ERROR: BKQ dict is missing in MatrixH')
                exit()

        if eig_opt:
            print('\nWavefunction optimization...')
            dic_bkq, quat = self.opt_eigenfunction_minimization(self, dic_bkq, self.basis.shape[0])
            print('...done')
            print('CFP rotation quaternion: ', quat)
        elif not eig_opt and cfp_angles is not None and len(cfp_angles)==3:
            dic_bkq = rota_LF(self.l, dic_bkq, *cfp_angles)
        elif not eig_opt and cfp_angles is not None and len(cfp_angles)>3:
            dict, coeff = read_DWigner_quat()
            dic_bkq = rota_LF_quat(self.l, dic_bkq, cfp_angles, dict, coeff)
        else:
            pass

        self.par_dict = {'F0':F0, 'F2':F2, 'F4':F4, 'F6':F6, 'zeta':zeta, 'k':k, 'dic_bkq':dic_bkq, 'field':field} 
        matrix = self.build_matrix(elem, F0, F2, F4, F6, zeta, k, dic_bkq, field, save_label, save_LF, save_matrix)

        w,v = diagonalisation(matrix)
        v = np.round(v, 16)  #numbers are saved as complex128 data type
        result = np.vstack((w,v))

        result = np.copy(result[:, result[0,:].argsort()])

        E = np.copy(result[0,:]).real

        print('Calculation result: \n' if wordy else "", end = "")
        print('E0: '+f'{min(w):.3f}'+'\n' if wordy else "", end = "")

        ####just pretty-print
        energy_print = np.around(E-min(E),8)
        energy_list, energy_count = np.unique(energy_print, return_counts=True)
        for i in range(len(energy_list)):
            if energy_count[i]!=1:
                deg_str = f' ({energy_count[i]})'
            else:
                deg_str = ''
            print(f' {energy_list[i]:10.3f}'+deg_str+'\n' if wordy else "", end = "")
        ####

        if ground_proj == True:
            print('\nGround state projection: ' if wordy else "", end = "")
            if self.ground==True:
                projected = projection_basis(result[1:,:], self.basis_l_JM, J_label=True)
                if wordy:
                    pprint(projected[1])  
            else:
                projected = projection_basis(result[1:,:], self.basis_l, J_label=True), projection_basis(result[1:,:], self.basis_l_JM, J_label=True)
                if wordy:
                    pprint(projected[0][1]) 
                    pprint(projected[1][1])

        if Orth==True:
            print('Orthogonality check... \n' if wordy else "", end = "")
            for i in range(len(w)):
                for j in range(len(w)):
                    if i != j:
                        check = np.abs(np.dot(np.conj(result[1:,i]).T, result[1:,j]))
                        if round(check, 5) != 0:
                            warnings.warn('Non-orthognal eigenvectors found')
                            print(i, j, check)
            print('...done\n' if wordy else "", end = "")
        if Norm==True:
            print('Normalization check... \n' if wordy else "", end = "")
            for i in range(len(w)):
                check = np.abs(np.dot(np.conj(result[1:,i]).T, result[1:,i]))
                if round(check, 5) != 1:
                    warnings.warn('Non-normalized eigenvectors found')
                    print(i, check)
            print('...done\n' if wordy else "", end = "")

        if return_proj:
            return result, projected
        else:
            return result

    # @cron
    def build_matrix(self, elem, F0=0, F2=0, F4=0, F6=0, zeta=0, k=1, dic_bkq = None, field = [0.,0.,0.], save_label=False, save_LF=False, save_matrix=False):
        """
        Builds the Hamiltonian matrix.

        Parameters:
        elem (list): The contributions to the Hamiltonian matrix among: 
                    'Hee'=interelectronic repulsion, 'Hso'=spin-orbit coupling, 
                    'Hcf'=crystal-field/ligand-field contribution, 'Hz'=Zeeman splitting.
        F0 (float, optional): The F0 Slater-Condon parameter in cm^{-1}. Default is 0.
        F2 (float, optional): The F2 Slater-Condon parameter in cm^{-1}. Default is 0.
        F4 (float, optional): The F4 Slater-Condon parameter in cm^{-1}. Default is 0.
        F6 (float, optional): The F6 Slater-Condon parameter in cm^{-1}. Default is 0.
        zeta (float, optional): The spin-orbit coupling parameter in cm^{-1}. Default is 0.
        k (float, optional): The orbital reduction factor. Default is 1.
        dic_bkq (dict, optional): A dictionary containing the crystal field parameters in Wybourne formalism. Default is None.
                                  dic_bkq = {'k':{'q':0.0, '-q':0.0}, ...}, with k=2,4,..2*l and q=-k,-k+1,...k-1,k.
        field (list, optional): The magnetic field vector in Tesla. Default is [0.,0.,0.].
        save_label (bool, optional): A flag indicating whether to save the labels of the states. Default is False.
        save_LF (bool, optional): A flag indicating whether to save the Hamiltonian matrix (without the Zeeman contribution). Default is False.
        save_matrix (bool, optional): A flag indicating whether to save the Hamiltonian matrix. Default is False.

        Returns:
        matrix (numpy.ndarray): The Hamiltonian matrix (complex).
        """

        F = F0, F2, F4, F6

        basis = self.basis
        dic_LS = self.dic_LS

        if self.closed==True:
            fac = -1.
        else:
            fac = 1.

        matrix = np.zeros((basis.shape[0],basis.shape[0]),dtype='complex128')
        if save_label:
            label_matrix = []
        if save_LF:
            LF_matrixs = np.zeros_like(matrix)

        for i in range(basis.shape[0]):
            statei = basis[i]
            Si = statei[0]/2.
            Li = statei[1]
            Ji = statei[2]/2.
            MJi = statei[3]/2.
            seni = statei[4]
            labeli = dic_LS[':'.join([f'{qq}' for qq in statei])]

            if save_label:
                label_matrix.append([Si*2,Li,Ji*2,MJi*2,seni])

            for j in range(0,i+1):
                statej = basis[j]
                Sj = statej[0]/2.
                Lj = statej[1]
                Jj = statej[2]/2.
                MJj = statej[3]/2.
                senj = statej[4]
                labelj = dic_LS[':'.join([f'{qq}' for qq in statej])]
                
                H = Hamiltonian([seni,Li,Si,senj,Lj,Sj,Ji,MJi,Jj,MJj], [labeli,labelj], self.conf, self.dic_cfp, self.tables, dic_LS, self.dic_LS_almost)  

                if 'Hee' in elem:
                    if Ji==Jj and MJi==MJj:
                        if Li == Lj and Si == Sj:
                            if self.l==3:
                                Hee = H.electrostatic_int(basis, *F, tab_ee = self.dic_ee)
                            else:
                                Hee = H.electrostatic_int(basis, *F)
                            matrix[i,j] += Hee
                            if save_LF:
                                LF_matrixs[i,j] += Hee
                                if i!=j:
                                    LF_matrixs[j,i] += np.conj(Hee)
                            if i != j:
                                matrix[j,i] += Hee

                if 'Hso' in elem:
                    if Ji==Jj and MJi==MJj:
                        Hso = fac*H.SO_coupling(zeta, k)
                        matrix[i,j] += Hso
                        if save_LF:
                            LF_matrixs[i,j] += Hso
                            if i!=j:
                                LF_matrixs[j,i] += np.conj(Hso)
                        if i != j:
                            matrix[j,i] += Hso

                if 'Hcf' in elem:
                    if Si==Sj:
                        Hcf = fac*H.LF_contribution(dic_bkq)
                        matrix[i,j] += Hcf
                        if save_LF:
                            LF_matrixs[i,j] += Hcf
                            if i!=j:
                                LF_matrixs[j,i] += np.conj(Hcf)
                        if i != j:
                            matrix[j,i] += np.conj(Hcf)

                if 'Hz' in elem:
                    if Li==Lj and Si==Sj and seni==senj:
                        Hz = H.Zeeman(field, k)
                        matrix[i,j] += Hz
                        if i != j:
                            matrix[j,i] += np.conj(Hz)

                else:
                    pass

        if save_label:
            np.savetxt('matrix_label.txt', np.array(label_matrix))
        if save_LF:
            np.save('matrix_LF', LF_matrixs, allow_pickle=True, fix_imports=False)
        if save_matrix:
            np.save('matrix', matrix, allow_pickle=True, fix_imports=False)

        return matrix

    @staticmethod
    def opt_eigenfunction_minimization(calc, dic_Bkq, shape):
        """
        Wavefunction optimization algorithm. It's based on the maximization of a single component of the ground state wavefunction
        using the dual anneling minimization algorithm implemented in scipy.

        Parameters:
        calc (object): The instance of the class.
        dic_Bkq (dict): A dictionary containing the crystal field parameters in Wybourne Default is None.
                        dic_bkq = {'k':{'q':0.0, '-q':0.0}, ...}, with k=2,4,..2*l and q=-k,-k+1,...k-1,k.
        shape (int): The shape of the ground state wavefunction.

        Returns:
        dic_rot (dict): A dictionary containing the rotated crystal field parameters.
        quat (list): A list containing the quaternion for the rotation of the crystal field parameters.
        """

        from scipy.optimize import dual_annealing
        from scipy.optimize import least_squares
        import time

        def use_nja_(calc, dic, wordy=False):
            result = calc.MatrixH(['Hcf'], **dic, eig_opt=False, wordy=wordy)
            ground_state = np.abs(result[1:,0])**2
            or_ground_state = np.sort(ground_state)
            return or_ground_state

        cycle = 0
        def target(quat, x):
            nonlocal cycle
            cycle += 1
            dic_rot = rota_LF_quat(3, dic_Bkq, quat, dict, coeff)
            dic = {'dic_bkq': dic_rot}
            ground_comp = use_nja_(calc, dic, wordy=False)*100
            print(f'Opt cycle: {cycle}, target: {np.sum((ground_comp-x)**2):.4e}, highest comp: {ground_comp[-1]:.2f}', end='\r')
            return np.sum((ground_comp-x)**2) 
        
        dict, coeff = read_DWigner_quat()

        x = np.zeros(shape)
        x[-1] = 100

        start_time = time.time()
        bounds = [(-1, 1), (-1, 1), (-1, 1), (-1, 1)]
        bounds_ls = ([-1, -1, -1, -1], [1, 1, 1, 1])
        result_global = dual_annealing(target, bounds, args=(x,), accept=-5.0, maxiter=500, maxfun=1000, seed=666, no_local_search=True)
        result = least_squares(target, result_global.x, bounds=bounds_ls, args=(x, ))
        end_time = time.time()
        elapsed_time = end_time - start_time
        print('\n'+f"Opt time: {elapsed_time:.2f} seconds"+'\n')

        dic_rot = rota_LF_quat(calc.l, dic_Bkq, result.x, dict, coeff)

        R = result.x
        quat = [R[-1], R[0], R[1], R[2]]

        return dic_rot, quat

    @staticmethod
    def opt_eigenfunction_grid(calc, dic_Bkq, shape):
        """
        Wavefunction optimization algorithm. It's based on a grid search, using REPULSION angles (evenly sampled on a sphere), 
        read from crystdat.py.

        Parameters:
        calc (object): The instance of the class.
        dic_Bkq (dict): A dictionary containing the crystal field parameters in Wybourne Default is None.
                        dic_bkq = {'k':{'q':0.0, '-q':0.0}, ...}, with k=2,4,..2*l and q=-k,-k+1,...k-1,k.
        shape (int): The shape of the ground state wavefunction.

        Returns:
        dic_rot (dict): A dictionary containing the rotated crystal field parameters.
        quat (list): A list containing the quaternion for the rotation of the crystal field parameters.
        """

        import time

        def use_nja_(calc, dic, wordy=False):
            result = calc.MatrixH(['Hcf'], **dic, eig_opt=False, wordy=wordy)
            ground_state = np.abs(result[1:,0])**2
            or_ground_state = np.sort(ground_state)
            return or_ground_state
        
        dict, coeff = read_DWigner_quat()

        x = np.zeros(shape)
        x[-1] = 1

        rep_cryst = np.array(crystdat.rep678_cryst)  #change here the name of the grid to select smaller or bigger grids from cristdat.py

        target_list = []

        start_time = time.time()

        for i in range(rep_cryst.shape[0]):

            angles = rep_cryst[i,:]
            a = angles[0]
            b = angles[1]

            r = scipy.spatial.transform.Rotation.from_euler('ZYZ', [a,b,0], degrees=True)
            R = r.as_quat()
            quat = [R[-1], R[0], R[1], R[2]]

            dic_rot = rota_LF_quat(calc.l, dic_Bkq, quat, dict, coeff)
            dic = {'dic_bkq': dic_rot}
            
            ground_comp = use_nja_(calc, dic, wordy=False)

            target = np.sum((np.abs(ground_comp-x))**2)
            target_list.append(target)

        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Opt time: {elapsed_time:.2f} seconds"+'\n')

        target_list = np.array(target_list)

        min_index = np.argmin(target_list)
        angles = rep_cryst[min_index,:]
        a = angles[0]
        b = angles[1]
        r = scipy.spatial.transform.Rotation.from_euler('ZYZ', [a,b,0], degrees=True)
        R = r.as_quat()
        quat = [R[-1], R[0], R[1], R[2]]
        dic_rot = rota_LF_quat(calc.l, dic_Bkq, quat, dict, coeff)

        return dic_rot, quat
