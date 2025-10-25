class Hamiltonian():
    """
    This class presents the routines for the computation of the Hamiltonian matrix element for the given integral. 

    References:
    1. R. Boca, "A handbook of magnetochemical formulae". Elsevier, 2012.
    2. R. Boca, "theoretical foundations of molecular magnetism" 1999. 
    3. E. Konig "Ligand Field Energy Diagrams" 2013.
    4. C. Goerller-Walrand, K. Binnemans, Handbook of Physics & Chemistry of Rare Earths, Vol 23, Ch 155, 1996.

    Attributes:
    v, L, S, v1, L1, S1 (float): The quantum numbers for the state.
    label1, label2 (str): The labels for the states.
    dic_LS (dict): A dictionary for the LS-coupling scheme.
    dic_LS_inv_almost (dict): An inverse dictionary for the LS-coupling scheme for the almost-closed-shell configuration.
    s (float): The spin quantum number.
    l (int): The one-electron orbital quantum number.
    n (int): The number of electrons in the almost-closed-shell configuration.
    N (int): The number of electrons in the configuration.
    cfp (CFP): A CFP object for the configuration.
    closed (bool): A flag indicating whether the configuration is almost-closed-shell or not.
    """

    def __init__(self, state, labels, conf, dic_cfp=None, tables=None, dic_LS=None, dic_LS_almost=None):
        #state = [v, L, S, v1, L1, S1, J, M, J1, M1]
        #labels = [label1, label2]
        """
        Initializes the Hamiltonian object with the given state, configuration, CFP dictionary, labels, and LS-coupling scheme dictionaries.

        The state is represented by a list of the form [v, L, S, v1, L1, S1, J, MJ, J1, MJ1], where 'v', 'L', 'S', 'J' and 'MJ' are the quantum numbers 
        for the initial state and 'v1', 'L1', 'S1', 'J1' and 'MJ1' are the quantum numbers for the final state. 
        The configuration is represented by a string of the form 'nl^x', where 'n' is the principal quantum number, 
        'l' is the azimuthal quantum number (represented as 'd' for l=2 and 'f' for l=3), 
        and 'x' is the number of electrons in the configuration.

        Parameters:
        state (list): The quantum numbers for the state.
        labels (list): The labels for the states.
        conf (str): The electron configuration.
        dic_cfp (dict): A dictionary containing the coefficients of fractional parentage (CFP) for the configuration.
        dic_LS (dict): A dictionary for the LS-coupling scheme.
        dic_LS_inv_almost (dict): An inverse dictionary for the LS-coupling scheme for the almost-closed-shell configuration.
        """

        self.v = state[0]
        self.L = state[1]
        self.S = state[2]
        self.v1 = state[3]
        self.L1 = state[4]
        self.S1 = state[5]
        self.J = state[6]
        self.M = state[7]
        self.J1 = state[8]
        self.M1 = state[9]
        self.label1 = labels[0]
        self.label2 = labels[1]

        self.conf = conf
        self.s = 0.5
        if conf[0]=='d':
            self.l = 2
        else:
            self.l = 3

        self.n = int(conf[1:])
        self.N = int(conf[1:])  #quelli veri
        self.closed = False
        if self.l == 2 and self.n>5:
            self.closed = True
            self.n = almost_closed_shells(conf)
        elif self.l == 2 and self.n<=5:
            pass
        elif self.l == 3 and self.n>7:
            self.closed = True
            self.n = almost_closed_shells(conf)
        elif self.l == 3 and self.n<=7:
            pass

        if dic_cfp is not None:
            self.TAB = False
            if self.closed==True:
                self.dic_LS_inv_almost = {}
                for k1,v1 in dic_LS_almost.items():
                    if v1 in self.dic_LS_inv_almost.keys():
                        self.dic_LS_inv_almost[v1] += [eval('['+k1.replace(':',',')+']')]
                    else:
                        self.dic_LS_inv_almost.update({v1:[eval('['+k1.replace(':',',')+']')]})
            else:
                self.dic_LS_inv_almost = None
            self.rme = RME(state, conf, dic_cfp, labels, dic_LS, self.dic_LS_inv_almost)
            self.dic_LS = dic_LS
        elif tables is not None and dic_cfp is None:
            self.TAB = True
            self.rme = tables
    
    def electrostatic_int(self, basis, F0=0, F2=0, F4=0, F6=0, evaluation=True, tab_ee=None):
        """
        Computes electron repulsion integrals <self.label1|Hee|self.label2> for the given basis set.
        Equations are taken from Boca 1999, "theoretical foundations of molecular magnetism" (Ch 8, p 518) (valid only for l2 conf)
        For the d^n configurations the equations are reported in Boca2012 (p 145 eq 4.66-4.69)

        Parameters:
        basis (numpy.ndarray): The basis set for the calculation.
        F0, F2, F4, F6 (float, optional): The F0, F2, F4, and F6 are the Slater-Condon parameters in cm^{-1}. F0 is irrelevant.
        evaluation (bool, optional): A flag indicating whether the parameters are symbolic or numerical. Default is True.
        tab_ee (dict, optional): A dictionary containing the electron-electron integrals. Default is None.

        Returns:
        integral (float): The calculated electron repulsion integral
        """

        def l_Ck_l1(l, k, l1):
            return (-1)**l*np.sqrt((2*l+1)*(2*l1+1))*Wigner_coeff.threej_symbol([[l, k, l1],[0, 0, 0]])

        def Vee(label1v, label1v1):
            #works only for d^n, f^1 and f^2

            idxv = term_label.index(label1v)
            Sv, Lv, vv = term_basis[idxv]
            Sv /= 2

            idxv1 = term_label.index(label1v1)
            Sv1, Lv1, vv1 = term_basis[idxv1]
            Sv1 /= 2

            ck_list = np.zeros(len(range(0,2*self.l+1,2)))
            for i,k in enumerate(range(0,2*self.l+1,2)):
                ck = 0
                if k!=0:
                    ck += 0.5*l_Ck_l1(self.l, k, self.l)**2
                    somma = 0
                    for ii, term in enumerate(term_basis):
                        if self.TAB == False:

                            label2 = term_label[ii]
                            S1, L1, v1 = term_basis[ii]  #2S, L, v
                            S1 /= 2

                            if label2[0]==label1v[0]:
                                somma += self.rme.Uk(k, n=self.N, v=vv, L=Lv, S=Sv, v1=v1, L1=L1, S1=S1, label1=label1v, label2=label2)*self.rme.Uk(k, n=self.N, v=vv1, L=Lv1, S=Sv1, v1=v1, L1=L1, S1=S1, label1=label1v1, label2=label2)  
                        else:

                            label2 = term_label[ii]
                            if label2[0]==label1v[0]:
                                try:
                                    somma += self.rme['U'+str(k)][label1v][label2]*self.rme['U'+str(k)][label1v1][label2]
                                except:
                                    somma += 0
                    somma *= 1/(2*Lv+1)
                    if self.v==self.v1:
                        somma -= self.n/(2*self.l+1)
                    ck *= somma
                else:
                    if self.v==self.v1:
                        ck = self.n*(self.n-1)/2
                if self.closed==True:
                    if self.v==self.v1:
                        ck -= (2*self.l+1-self.n)/(2*self.l+1)*l_Ck_l1(self.l, k, self.l)**2
                else:
                    pass
                ck_list[i] = ck

            return ck_list
           
        #=======================================================================

        if evaluation==False:
            F0, F2, F4, F6 = sympy.symbols("F0, F2, F4, F6")
            coeff = [F0, F2, F4, F6]
        else:
            coeff = [F0, F2, F4, F6]

        if tab_ee is not None:
            coeff_ee = tab_ee.copy()

        term_basis = np.array(terms_basis(self.conf[0]+str(self.n)))
        term_label = terms_labels(self.conf[0]+str(self.n))
        ck_coeff = np.zeros(len(range(0,2*self.l+1,2)))
        for i, term in enumerate(term_basis):
            if tab_ee is not None:
                ck_coeff += coeff_ee[term_label[i]][term_label[i]]*(term[0]+1)*(2*term[1]+1)
            else:
                ck_coeff += np.array(Vee(term_label[i], term_label[i]), dtype='float64')*(term[0]+1)*(2*term[1]+1)

        ck_coeff /= basis.shape[0]

        integral = 0
        label1v = self.label1
        label1v1 = self.label2
        if tab_ee is not None:
            try:
                ck = np.array(coeff_ee[label1v][label1v1])
            except:
                ck = np.zeros(4)
        else:
            ck = Vee(label1v, label1v1)
        for i,_ in enumerate(range(0,2*self.l+1,2)):
            if label1v==label1v1:
                integral += (ck[i] - ck_coeff[i])*coeff[i]
            else:
                integral += ck[i]*coeff[i]
        
        return integral

    def SO_coupling(self, zeta, k=1, evaluation=True):
        """
        Computes SOC integrals <self.label1|Hso|self.label2>.
        Equations are taken from Konig & Kremer (Ch 2, p 11, eq 2.82)

        Parameters:
        zeta (float): The SOC parameter in cm^{-1}.
        k (int, optional): The orbital reduction factor. Default is 1.
        evaluation (bool, optional): A flag indicating whether the parameters are symbolic or numerical. Default is True.

        Returns:
        integral (float): The calculated SOC integral.
        """

        if evaluation==False:
            zeta, k = sympy.symbols("zeta, k")
        else:
            pass

        pref = zeta*k*(-1)**(self.J+self.L+self.S1)*np.sqrt(self.l*(self.l+1)*(2*self.l+1))
        coeff = Wigner_coeff.sixj_symbol([[self.L, self.L1,1],[self.S1,self.S, self.J]])
        if self.TAB == False:
            rme_V1k = self.rme.V1k()
        else:
            try:
                rme_V1k = self.rme['V11'][self.label1][self.label2]
            except:
                rme_V1k = 0

        return pref*coeff*rme_V1k

    def LF_contribution(self, dic_kq, evaluation=True):
        """
        Computes the LF/CF contribution to the Hamiltonian matrix element <self.label1|HCF/LF|self.label2> using Bkq in Wybourne formalism.
        Equations are taken from C. Goerller-Walrand, K. Binnemans, Handbook of Physics & Chemistry of Rare Earths, Vol 23, Ch 155, (1996)
        or eq 3 of Rudowicz and Karbowiak Coordination chemistry reviews 287 (2015) 28-63.

        Parameters:
        dic_kq (dict): A dictionary containing the Bkq coefficients in cm^{-1}, e.g. dic_bkq = {'2':{'-1':0.0}} where k=2 and q=-1. 
        The complete (complex) coefficient corresponds to Bkq + iBk-q. 
        evaluation (bool, optional): A flag indicating whether the parameters are symbolic or numerical. Default is True.

        Returns:
        result (complex): The calculated LF/CF contribution to the Hamiltonian matrix element.
        """

        if evaluation==False:
            for k in range(2,2*self.l+1,2):
                for q in range(-k,k+1,1):
                    dic_kq[f'{k}'][f'{q}'] = sympy.Symbol(f"B{k}{q}")

        else:
            for k in range(2,2*self.l+1,2):
                for q in range(-k,k+1,1):
                    if f'{k}' in dic_kq.keys() and f'{q}' in dic_kq[f'{k}'].keys():
                        pass
                    else:
                        dic_kq[f'{k}'][f'{q}'] = 0.0

        result = 0
        for k in range(2,2*self.l+1,2):
            ckk = (-1)**self.l*(2*self.l+1)*Wigner_coeff.threej_symbol([[self.l, k, self.l],[0,0,0]])
            if self.TAB == False:
                Uk = self.rme.Uk(k)
            else:
                try:
                    Uk = self.rme['U'+str(k)][self.label1][self.label2]
                except:
                    Uk = 0
           
            coeff2 = Wigner_coeff.sixj_symbol([[self.J, self.J1, k],[self.L1, self.L, self.S]])
           
            integral = 0
            for q in range(0,k+1):
                if q==0:
                    coeff0 = Wigner_coeff.threej_symbol([[self.J, k, self.J1],[-self.M,0,self.M1]])
                    Yk0 = (-1)**(2*self.J + self.L1 + self.S - self.M + k)*np.sqrt((2*self.J+1)*(2*self.J1+1))*ckk*coeff0*coeff2*Uk
                    integral += dic_kq[f'{k}'][f'{0}']*Yk0
                else:
                    coeffp = Wigner_coeff.threej_symbol([[self.J, k, self.J1],[-self.M,q,self.M1]])
                    coeffm = Wigner_coeff.threej_symbol([[self.J, k, self.J1],[-self.M,-q,self.M1]])
                    Ckqp = (-1)**(2*self.J + self.L1 + self.S - self.M + k)*np.sqrt((2*self.J+1)*(2*self.J1+1))*ckk*coeffp*coeff2*Uk
                    Ckqm = (-1)**(2*self.J + self.L1 + self.S - self.M + k)*np.sqrt((2*self.J+1)*(2*self.J1+1))*ckk*coeffm*coeff2*Uk
                    if eval:
                        integral += dic_kq[f'{k}'][f'{q}']*(Ckqm+(-1)**q*Ckqp) + 1j*dic_kq[f'{k}'][f'{-q}']*(Ckqm - (-1)**q*Ckqp)
                    else:
                        integral += dic_kq[f'{k}'][f'{q}']*(Ckqm+(-1)**q*Ckqp) + sympy.I*dic_kq[f'{k}'][f'{-q}']*(Ckqm - (-1)**q*Ckqp)

            result += integral

        return result

    def Zeeman(self, field=np.array([0.,0.,0.]), k=1, evaluation=True, MM=False):
        """
        Computes the Zeeman contribution to the Hamiltonian matrix element <self.label1|HZeeman|self.label2> using the magnetic field.
        Equations are taken from Boca2012, "A handbook of magnetochemical formulae" (p 588).

        Parameters:
        field (numpy.ndarray): The magnetic field vector in T.
        k (int, optional): The rank of the tensor operator. Default is 1.
        evaluation (bool, optional): A flag indicating whether the parameters are symbolic or numerical. Default is True.
        MM (bool, optional): A flag indicating whether to return the L1 and S1 contributions. Default is False.

        Returns:
        int (complex): The calculated Zeeman contribution to the Hamiltonian matrix element.
        """

        if evaluation ==False:
            Bx, By, Bz = sympy.symbols("Bx, By, Bz")
        else:
            Bx, By, Bz = np.array(field)

        Bohr = 0.4668604  #conv from BM to cm-1
        ge = 2.0023

        Bq = [-np.sqrt(0.5)*(Bx+1j*By), Bz, np.sqrt(0.5)*(Bx-1j*By)]  # +1, 0, -1
        L1 = []
        S1 = []

        pre = np.sqrt((2*self.J+1)*(2*self.J1+1))

        L1q = k*np.sqrt(self.L*(self.L+1)*(2*self.L+1))*(-1)**(self.L+self.S+self.J+1)*Wigner_coeff.sixj_symbol([[self.L, self.J, self.S],[self.J1, self.L, 1]])#Wigner_coeff.sixj_symbol([[self.J, 1, self.J1],[self.L, self.S, self.L]])#Wigner_coeff.sixj_symbol([[self.J1, self.J, 1],[self.L, self.L, self.S]])

        S1q = ge*np.sqrt(self.S*(self.S+1)*(2*self.S+1))*(-1)**(self.L+self.L1+self.S1*2+self.J+self.J1)*(-1)**(self.L+self.S+self.J+1)*Wigner_coeff.sixj_symbol([[self.S, self.J, self.L],[self.J1, self.S, 1]])#Wigner_coeff.sixj_symbol([[self.J, 1, self.J1],[self.S, self.L, self.S]])#Wigner_coeff.sixj_symbol([[self.J1, self.J, 1],[self.S, self.S, self.L]])

        rme = pre*(L1q + S1q)

        int = 0
        for i,q in enumerate(range(-1,2,1)):
            preq = (-1)**(self.J-self.M)*Wigner_coeff.threej_symbol([[self.J, 1, self.J1],[-self.M, q, self.M1]])
            int += (-1)**q*Bq[i]*preq*rme*Bohr

            L1.append(preq*pre*L1q)
            S1.append(preq*pre*S1q)

        if not MM:
            return int
        else:
            return L1, S1
