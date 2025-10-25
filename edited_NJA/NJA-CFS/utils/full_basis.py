import numpy as np

def Full_basis(conf):
    """
    Produces the complete basis set for a given configuration, saved in different format:
    1. basis --> array: n. microstates x [2S, L, 2J, 2M, sen (,count)] 
    2. dic_LS --> dict: '2S:L:2J:2M:sen:count': label as N. and K.
    3. basis_l --> list: 2S+1 L (J)
    4. basis_l_JM --> list: 2S+1 L (J) MJ

    Parameters:
    conf (str): The electron configuration.

    Returns:
    basis (numpy.ndarray): The complete basis set for the configuration.
    dic_LS (dict): A dictionary for the LS-coupling scheme.
    basis_l (numpy.ndarray): The labels for the states.
    basis_l_JM (numpy.ndarray): The labels for the states (with MJ).
    """

    n_el = int(conf[1:])

    if conf[0]=='d' and n_el>5:
        n_el = 10-n_el
    elif conf[0]=='f' and n_el>7:
        n_el = 14-n_el
    else:
        pass

    if conf[0]=='d':
        n_freeion_SL = [1,5,8,16,16]
    else:
        n_freeion_SL = [1,7,17,47,73,119,119]

    basis = []
    basis_l = []
    basis_l_JM = []
    dic_LS = {}
    count=0

    for st in range(0,n_freeion_SL[n_el-1]):  #run over free ion states
        term_base = terms_basis(conf)[st]
        TwoS = term_base[0]
        L = term_base[1]
        sen = term_base[2]
        label_l = terms_labels(conf)[st]

        if len(label_l)>2:
            sen_str = label_l[-1]
        else:
            sen_str = ''
        if sen_str=='10':
            sen_str = '0'

        J1 = np.abs(2*L-TwoS)
        J2 = 2*L + TwoS

        for TwoJ in range(J1,J2+2,2):
            if TwoJ%2==0:
                J_str = str(int(TwoJ/2))
            else:
                J_str = str(int(TwoJ))+'/2'
            for TwoMJ in range(-TwoJ,TwoJ+2,2):
                if conf[0]=='f':
                    if sen_str=='':
                        count = 0
                    else:
                        if sen_str=='0':
                            count = 10
                        else:
                            count = int(sen_str)
                    lista = [TwoS, L, TwoJ, TwoMJ, sen, count]
                else:
                    lista = [TwoS, L, TwoJ, TwoMJ, sen]
                basis.append(lista)
                dic_LS[':'.join(f'{n}' for n in lista)] = label_l
                if TwoMJ%2==0:
                    MJ_str = str(int(TwoMJ/2))
                else:
                    MJ_str = str(int(TwoMJ))+'/2'
                basis_l.append(str(TwoS+1)+state_legend(str(L), inv=True)+sen_str+' ('+J_str+')')  
                basis_l_JM.append(str(TwoS+1)+state_legend(str(L), inv=True)+sen_str+' ('+J_str+') '+MJ_str)

    basis = np.array(basis)
    basis_l = np.array(basis_l)
    basis_l_JM = np.array(basis_l_JM)

    return basis, dic_LS, basis_l, basis_l_JM
