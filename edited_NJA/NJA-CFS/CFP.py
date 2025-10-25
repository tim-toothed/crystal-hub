import numpy as np

class CFP():
    """
    Handles the coefficients of fractional parentage (cfp) for a given electron configuration.

    References:
    1. Nielson - Koster in "Spectroscopic Coefficients for the p^n, d^n and f^n Configurations" (1963)

    Attributes:
    dic_LS_inv_almost (dict): An inverse dictionary for the LS-coupling scheme for the almost-closed-shell configuration.
    l (int): The one-electron orbital quantum number.
    n (int): The number of electrons in the almost-closed-shell configuration.
    N (int): The number of electrons in the configuration.
    closed (bool): A flag indicating whether the configuration is almost-closed-shell or not.
    """

    def __init__(self, conf, dic_cfp, dic_LS_inv_almost=None):
        """
        Initializes the CFP object with the given configuration, CFP dictionary, and optional LS-coupling scheme dictionary.

        The configuration is represented by a string of the form 'l^x', where 'l' is orbital momentum quantum number 
        (represented as 'd' for l=2 and 'f' for l=3), and 'x' is the number of electrons in the configuration.

        Parameters:
        conf (str): The electron configuration.
        dic_cfp (dict): A dictionary containing the coefficients of fractional parentage (CFP) for the configuration.
        dic_LS_inv_almost (dict, optional): An inverse dictionary for the LS-coupling scheme for the almost-closed-shell configuration. Default is None.
        """

        if conf[0]=='d':
            self.l = 2
        else:
            self.l = 3

        self.n = int(conf[1:])  
        self.N = int(conf[1:])  

        self.dic_cfp = dic_cfp

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

        self.dic_LS_inv_almost = dic_LS_inv_almost


    def cfp(self, v, L, S, name):
        #The equation for the almost-closed-shell parameters calculation is given in 
        # Nielson - Koster in "Spectroscopic Coefficients for the p^n, d^n and f^n Configurations" (1963)
        """
        Returns the cfps for a given state, identified by its L, S quantum number and the name of the state.

        Parameters:
        L (int): Orbital angular momentum quantum number
        S (float): Spin angular momentum quantum number
        name (str): The name of the state (from terms_labels()).

        Returns:
        cfp_list (numpy.ndarray): The list of CFPs for the given state.
        """

        dic = self.dic_cfp 

        if self.closed==True:  
            cfp_list = []
            for keys in dic.keys():
                values_list = [[key, float(val)] for key, val in dic[keys].items()]
                values = sum(values_list, [])
                for i in range(0,len(values),2):
                    zz=1
                    if self.N == 2*self.l:
                        zz = (-1)**(v-1/2)
                    if values[i]==name:
                        term = self.dic_LS_inv_almost[keys][0] 
                        Sk = term[0]/2.
                        Lk = term[1]
                        N = self.N-1
                        cfp_value = values[i+1]/(zz*(-1)**(Sk+S+Lk+L-self.l-0.5)*np.sqrt((N+1)*(2*S+1)*(2*L+1)/((4*self.l+2-N)*(2*Sk+1)*(2*Lk+1))))  
                        cfp_list.append([keys,cfp_value])
        else:
            cfp_list = [[key, float(val)] for key, val in dic[name].items()]

        cfp_list = np.array(cfp_list, dtype='object')

        return cfp_list
