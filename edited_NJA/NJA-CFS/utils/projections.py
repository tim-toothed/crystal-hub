import numpy as np
import warnings

from .tables_legends_conventions import terms_labels

#########PROJECTIONS###########

def projection(basis2, labels, basis1, energy, bin=1e-4, min_contr=10):
    #just implementation of projection operator

    matrix_coeff = np.zeros_like(basis1)
    matrix_coeff_re2 = np.zeros_like(basis1, dtype='float64')
    states = []
    states_red = {}

    for i in range(basis2.shape[0]):
        states.append([])
        states_red[i+1] = {}
        for j in range(basis1.shape[0]):
            matrix_coeff[j,i] = np.dot(basis1[:,j].T,basis2[:,i])
            matrix_coeff_re2[j,i] = np.abs(np.dot(np.conj(basis1[:,j]).T,basis2[:,i]))**2
            stato = [f'{round(energy[j],3)}', matrix_coeff_re2[j,i]]
            states[i].append(stato)
            key, value = stato[0], stato[1]
            if value>bin:
                if key in states_red[i+1].keys():
                    states_red[i+1][key] += value
                else:
                    states_red[i+1][key] = value
            else:
                pass
        # print(matrix_coeff[:,i])   #combinazione lineare per basis2[:,i]
        # print(matrix_coeff_re2[:,i])
        # exit()
        tot = sum(states_red[i+1].values())
        if round(tot,2) != 1:
            warnings.warn('The expantion coefficient do not sum to 1')
            print(tot)
        for key, value in states_red[i+1].items():
            states_red[i+1][key] = value*100/tot
        sortato = sorted(states_red[i+1].items(), key=lambda x:x[1], reverse=True)  #sort dict on values
        states_red[i+1] = dict(sortato)

    # print(states_red)

    states_red_2 = {}   #in questo ci vanno solo i livelli con percentuali superiori o uguali a min_contr
    for key1 in states_red.keys():
        states_red_2[key1] = {}
        for key2, value in states_red[key1].items():
            if value>=min_contr:
                states_red_2[key1][key2] = value

    # print(states_red_2)

    return states_red_2

def projection_basis(basis2, labels, bin=1e-5, J_label=False):
    #just implementation of projection operator

    basis1 = np.eye(basis2.shape[0])
    matrix_coeff = np.zeros_like(basis1, dtype='complex128')
    matrix_coeff_re2 = np.zeros_like(basis1, dtype='float64')
    states = []
    states_red = {}
    for i in range(basis2.shape[0]):
        states.append([])
        states_red[i+1] = {}
        for j in range(basis1.shape[0]):
            matrix_coeff[j,i] = np.dot(basis1[:,j].T,basis2[:,i])
            matrix_coeff_re2[j,i] = np.abs(np.dot(np.conj(basis1[:,j]).T,basis2[:,i]))**2
            stato = [labels[j], matrix_coeff_re2[j,i]]
            states[i].append(stato)
            if J_label==True:
                key, value = stato[0], stato[1]
            else:
                key, value = stato[0][:stato[0].index(' (')], stato[1]
            if value>bin:
                if key in states_red[i+1].keys():
                    states_red[i+1][key] += value
                else:
                    states_red[i+1][key] = value
            else:
                pass
        tot = sum(states_red[i+1].values())
        if round(tot,2) != 1:
            warnings.warn('The expantion coefficient do not sum to 1')
            print(tot)
        for key, value in states_red[i+1].items():
            states_red[i+1][key] = value*100/tot
        sortato = sorted(states_red[i+1].items(), key=lambda x:x[1], reverse=True)  #sort dict on values
        states_red[i+1] = dict(sortato)

    return states_red

def the_highest_MJ(proj, state_list):

    cont = np.zeros_like(state_list)
    for key,value in proj.items():
        for i in range(len(cont)):
            end = key.find(')')
            M_key = np.abs(eval(key[end+1:]))
            if state_list[i]==M_key:
                cont[i] += value

    M = state_list[cont.argmax()]
    perc = max(cont)

    return M, perc

def the_highest_L(proj, conf):
    
    state_list = terms_labels(conf)
    
    cont = np.zeros(len(state_list))
    for key,value in proj.items():
        
        for i in range(len(cont)):
            try:
                L_key = key[:(key.index('(')-1)]
            except:
                L_key = key
            if state_list[i]==L_key:
                cont[i] += value

    L = state_list[cont.argmax()]
    perc = max(cont)

    return L, perc
