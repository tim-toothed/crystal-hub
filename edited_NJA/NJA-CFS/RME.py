class RME():  #RME(CFP)
    #from E. Konig "Ligand Field Energy Diagrams" (eq 2.85,2.87)
    #or from Boca, "theoretical foundations of molecular magnetism" (Ch 8, p 516)
    """
    This class is used to perform a RME calculation for a given electron configuration. 
    The electron configuration is represented by a string of the form 'l^x', where 'l' is the azimuthal quantum number 
    (represented as 'd' for l=2 and 'f' for l=3), and 'x' is the number of electrons in the configuration.

    References:
    1. E. Konig "Ligand Field Energy Diagrams" (eq 2.85,2.87)
    2. R. Boca, "theoretical foundations of molecular magnetism" (Ch 8, p 516)

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

    def __init__(self, state, conf, dic_cfp, labels, dic_LS, dic_LS_inv_almost):
        #state = [v, L, S, v1, L1, S1]
        #WARNING: the dic_LS is the one from the corresponding configuration in case of almost closed-shell + 1
        """
        Initializes the RME object with the given state, configuration, CFP dictionary, labels, and LS-coupling scheme dictionaries.

        The state is represented by a list of the form [v, L, S, v1, L1, S1], where 'v', 'L', and 'S' are the quantum numbers for the initial state and 'v1', 'L1', and 'S1' are the quantum numbers for the final state. The configuration is represented by a string of the form 'nl^x', where 'n' is the principal quantum number, 'l' is the azimuthal quantum number (represented as 'd' for l=2 and 'f' for l=3), and 'x' is the number of electrons in the configuration.

        Parameters:
        state (list): The quantum numbers for the state.
        conf (str): The electron configuration.
        dic_cfp (dict): A dictionary containing the coefficients of fractional parentage (CFP) for the configuration.
        labels (list): The labels for the states.
        dic_LS (dict): A dictionary for the LS-coupling scheme.
        dic_LS_inv_almost (dict): An inverse dictionary for the LS-coupling scheme for the almost-closed-shell configuration.
        """

        self.v = state[0]
        self.L = state[1]
        self.S = state[2]
        self.v1 = state[3]
        self.L1 = state[4]
        self.S1 = state[5]
        self.label1 = labels[0]
        self.label2 = labels[1]
        self.dic_LS = dic_LS
        self.dic_LS_inv_almost = dic_LS_inv_almost

        #super().__init__(self, conf)
        self.s = 0.5
        if conf[0]=='d':
            self.l = 2
        else:
            self.l = 3

        self.n = int(conf[1:])  #electrons closed conf
        self.N = int(conf[1:])  #n. of electrons true
        self.cfp = CFP(conf, dic_cfp, dic_LS_inv_almost) #cfp corrected for almost-closed-shell

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


    def Uk(self, k, n=None, v=None, L=None, S=None, v1=None, L1=None, S1=None, label1=None, label2=None):   
        """
        Calculate the Uk coefficient for a given set of quantum numbers and labels.

        This method calculates the Uk coefficient, which is used in the calculation of reduced matrix elements (RMEs) in atomic physics. 
        The quantum numbers and labels for the initial and final states are either provided as arguments or taken from the RME object itself. 
        The calculation involves summing over the coefficients of fractional parentage (CFP) for the initial and final states 
        and the 6-j symbol of the associated angular momenta.

        Parameters:
        k (int): The rank of the tensor operator.

        Returns:
        Uk (float): The calculated Uk coefficient.
        """

        if n is None:
            n = self.N
        if v is None or L is None or S is None or label1 is None:
            v, L, S, label1 = self.v, self.L, self.S, self.label1
        if v1 is None or L1 is None or S1 is None or label2 is None:
            v1, L1, S1, label2 = self.v1, self.L1, self.S1, self.label2

        pref = n*(2*L+1)**0.5*(2*L1+1)**0.5

        cfp_list = self.cfp.cfp(v, L, S, label1)
        cfp_list1 = self.cfp.cfp(v1, L1, S1, label2)

        somma = 0
        for i in range(len(cfp_list)):
            for j in range(len(cfp_list1)):
                if cfp_list[i,0] == cfp_list1[j,0]:
                    L_parent = state_legend(cfp_list[i,0][1])
                    matrix = [[L, L1, k],[self.l, self.l, L_parent]]
                    somma +=  cfp_list[i,-1]*cfp_list1[j,-1]*(-1)**(L_parent+L+self.l+k)*Wigner_coeff.sixj_symbol(matrix)

        if self.closed == True:
            return (-(-1)**k)*pref*somma
        else:
            return pref*somma


    def V1k(self, k=1): #[h^2]
        """
        Calculate the V1k coefficient for a given rank of the tensor operator.

        This method calculates the V1k coefficient, which is used in the calculation of reduced matrix elements (RMEs) in atomic physics. The rank of the tensor operator is provided as an argument. The calculation involves summing over the coefficients of fractional parentage (CFP) for the initial and final states and the 6-j symbols of the associated angular momenta.

        Parameters:
        k (int, optional): The rank of the tensor operator. Default is 1.

        Returns:
        V1k (float): The calculated V1k coefficient.
        """

        pref = self.N*((self.s*(self.s+1)*(2*self.s+1)))**0.5*((2*self.S+1)*(2*self.L+1)*(2*self.S1+1)*(2*self.L1+1))**0.5

        cfp_list = self.cfp.cfp(self.v, self.L, self.S, self.label1)  #cfps are already corrected for almost_closed_shells
        cfp_list1 = self.cfp.cfp(self.v1, self.L1, self.S1, self.label2)

        somma = 0
        for i in range(len(cfp_list)):
            for j in range(len(cfp_list1)):
                if cfp_list[i,0] == cfp_list1[j,0]:
                    L_parent = state_legend(cfp_list[i,0][1])
                    S_parent = (int(cfp_list[i,0][0])-1)/2
                    matrix1 = [[self.S, self.S1, 1],[self.s, self.s, S_parent]]
                    matrix2 = [[self.L, self.L1, k],[self.l, self.l, L_parent]]
                    somma += cfp_list[i,-1]*cfp_list1[j,-1]*(-1)**(L_parent+S_parent+self.S+self.L+self.s+self.l+k+1)*Wigner_coeff.sixj_symbol(matrix1)*Wigner_coeff.sixj_symbol(matrix2)

        if self.closed == True:
            return (-1)**k*pref*somma
        else:
            return pref*somma
