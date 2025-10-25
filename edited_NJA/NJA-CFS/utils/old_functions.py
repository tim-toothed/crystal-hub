import numpy as np
import warnings

#OTHER (old functions, not much used)

def from_Vint_to_Bkq(dic_V, conf):
    #conversion V matrix element to the Ckq coefficents of the ligand field expanded as spherical harmonics
    #see Gerloch & McMeeking (1975) (Table 2)
    #or OctoYot f_e_LF.f90 subroutine AOMmatrixD()
    #elements order: z2, yz, xz, xy, x2-y2

    if conf[0]=='d':
        l = 2
        dic_ckq = {
                  '0':{'0':2./5.*np.sqrt(np.pi)*(dic_V['11']+dic_V['22']+dic_V['33']+dic_V['44']+dic_V['55'])},
                  '2':{'0':np.sqrt(np.pi/5.)*(2.*dic_V['11'] + dic_V['22'] + dic_V['33'] - 2.*dic_V['44'] - 2.*dic_V['55']),
                       '1':-np.sqrt(4.*np.pi/5.)*( np.sqrt(3)/np.sqrt(2)*(dic_V['42'] + dic_V['53']) + dic_V['31']/np.sqrt(2)),
                       '-1':np.sqrt(4.*np.pi/5.)*( -np.sqrt(3)/np.sqrt(2)*(dic_V['52'] - dic_V['43']) + dic_V['21']/np.sqrt(2)),
                       '2':-np.sqrt(4.*np.pi/5.)*(np.sqrt(2)*dic_V['51'] + np.sqrt(3)/(2.*np.sqrt(2))*(dic_V['22']-dic_V['33'])),
                       '-2':-np.sqrt(4.*np.pi/5.)*(-np.sqrt(2)*dic_V['41'] + np.sqrt(3)/np.sqrt(2)*dic_V['32'])},
                   '4':{'0':np.sqrt(np.pi)/5.*(6.*dic_V['11']-4.*dic_V['22']-4.*dic_V['33']+dic_V['44']+dic_V['55']),
                        '1':2*np.sqrt(2.*np.pi/5.)*(-np.sqrt(3)/np.sqrt(2)*dic_V['31']+1./(2.*np.sqrt(2))*(dic_V['42'] + dic_V['53'])),
                        '-1':2*np.sqrt(2.*np.pi/5.)*(np.sqrt(3)/np.sqrt(2)*dic_V['21']+1./(2.*np.sqrt(2))*(dic_V['52'] - dic_V['43'])),
                        '2':2.*np.sqrt(2.*np.pi/5.)*(np.sqrt(3)/2.*dic_V['51']+(dic_V['33']-dic_V['22'])/2.),
                        '-2':2.*np.sqrt(2.*np.pi/5.)*(-np.sqrt(3)/2.*dic_V['41']-dic_V['32']),
                        '3':np.sqrt(7*np.pi/5.)*(dic_V['42']-dic_V['53']),
                        '-3':np.sqrt(7*np.pi/5.)*(dic_V['43']+dic_V['52']),
                        '4':np.sqrt(7.*np.pi/10.)*(dic_V['55']-dic_V['44']),
                        '-4':-np.sqrt(7.*np.pi/5.)*np.sqrt(2)*dic_V['54']}
                    }
        
    #conversion V matrix element to the Ckq coefficents of the ligand field expanded as spherical harmonics
    #given by Urland, Chem.Phys.14, 393,(1976). Table 3
    #or OctoYot f_e_LF.f90 subroutine AOMmatrixF()
    #                                        elements order: |sigma>, |piS>, |piC>, |deltaS>, |deltaC>, |phiS>, |phiC>
    #                                             (in ORCA):     0      -1     1       -2        2        -3      3
    #(For a definition of these orbitals: see Harnung & Schaffer, Struct&Bond,12,201,(1972))


    elif conf[0]=='f':
        l = 3
        dic_ckq = {
                  '0':{'0':(2./7.)*np.sqrt(np.pi)*(dic_V['11']+dic_V['22']+dic_V['33']+dic_V['44']+dic_V['55']+dic_V['66']+dic_V['77'])},
                  '2':{'0':(2./7.)*np.sqrt(5*np.pi)*dic_V['11'] + (3./14.)*np.sqrt(5*np.pi)*(dic_V['22']+dic_V['33']) - (5./14.)*np.sqrt(5*np.pi)*(dic_V['66']+dic_V['77']),
                       '1': - (1./7.)*np.sqrt(5*np.pi)*dic_V['31'] + (5./14.)*np.sqrt(3*np.pi)*(-dic_V['42']-dic_V['53']) + (5./14.)*np.sqrt(5*np.pi)*(-dic_V['64']-dic_V['75']),
                       '-1': (1./7.)*np.sqrt(5*np.pi)*dic_V['21'] + (5./14.)*np.sqrt(3*np.pi)*(dic_V['43']-dic_V['52']) + (5./14.)*np.sqrt(5*np.pi)*(dic_V['65']-dic_V['74']),
                       '2': (1./7.)*np.sqrt(30*np.pi)*(-dic_V['22']/2 + dic_V['33']/2) + (5./7.)*np.sqrt(2*np.pi)*(-dic_V['51']) + (5./7.)*np.sqrt(np.pi/2)*(-dic_V['62']-dic_V['73']),
                       '-2': (1./7.)*np.sqrt(30*np.pi)*(-dic_V['32']) + (5./7.)*np.sqrt(2*np.pi)*(dic_V['41']) + (5./7.)*np.sqrt(np.pi/2)*(dic_V['63']-dic_V['72'])},
                   '4':{'0': -np.sqrt(np.pi)*(dic_V['44']+dic_V['55']) + (1./7.)*np.sqrt(np.pi)*(6*dic_V['11'] + dic_V['22'] + dic_V['33'] + 3*dic_V['66'] + 3*dic_V['77']),
                        '1': (1./7.)*np.sqrt(30*np.pi)*(-dic_V['31']) + (4./7.)*np.sqrt(2*np.pi)*(-dic_V['42']-dic_V['53']) + (1./7.)*np.sqrt(30*np.pi)*(dic_V['64']+dic_V['75']),
                        '-1': (1./7.)*np.sqrt(30*np.pi)*(dic_V['21']) + (4./7.)*np.sqrt(2*np.pi)*(dic_V['43']-dic_V['52']) + (1./7.)*np.sqrt(30*np.pi)*(-dic_V['65']+dic_V['74']),
                        '2': (2./7.)*np.sqrt(10*np.pi)*(-dic_V['22']/2 + dic_V['33']/2) + (1./7.)*np.sqrt(6*np.pi)*(-dic_V['51']) + (3./7.)*np.sqrt(6*np.pi)*(dic_V['62']+dic_V['73']),
                        '-2': (2./7.)*np.sqrt(10*np.pi)*(-dic_V['32']) + (1./7.)*np.sqrt(6*np.pi)*(dic_V['41'] + 3*(-dic_V['63']+dic_V['72'])),
                        '3': np.sqrt((2./7.)*np.pi)*(dic_V['42']-dic_V['53']) + 3*np.sqrt((2./7.)*np.pi)*dic_V['71'],
                        '-3': np.sqrt((2./7.)*np.pi)*(dic_V['43']+dic_V['52']) - 3*np.sqrt((2./7.)*np.pi)*dic_V['61'],
                        '4': np.sqrt((10./7.)*np.pi)*(-dic_V['44']/2 + dic_V['55']/2) + np.sqrt((6./7.)*np.pi)*(dic_V['62']-dic_V['73']),
                        '-4': np.sqrt((10./7.)*np.pi)*(-dic_V['54']) + np.sqrt((6./7.)*np.pi)*(dic_V['63']+dic_V['72'])},
                   '6':{'0':(1./7.)*np.sqrt(13*np.pi)*(2*dic_V['11'] - (3./2.)*(dic_V['22']+dic_V['33']) + (3./5.)*(dic_V['44']+dic_V['55']) - (1./10.)*(dic_V['66']+dic_V['77'])),
                        '1':np.sqrt((13./7.)*np.pi)*(-dic_V['31'] + (1./2.)*np.sqrt(3./5.)*(dic_V['42']+dic_V['53']) + (1./10.)*(-dic_V['64']-dic_V['75'])),
                        '-1':np.sqrt((13./7.)*np.pi)*(dic_V['21'] + (1./2.)*np.sqrt(3./5.)*(-dic_V['43']+dic_V['52']) + (1./10.)*(dic_V['65']-dic_V['74'])),
                        '2':np.sqrt((13./7.)*np.pi)*(-(1./2.)*np.sqrt(3./5.)*(dic_V['22']-dic_V['33']) + (4./5.)*dic_V['51'] - (1./5.)*(dic_V['62']+dic_V['73'])),
                        '-2':np.sqrt((13./7.)*np.pi)*(-np.sqrt(3./5.)*dic_V['32'] - (4./5.)*dic_V['41'] + (1./5.)*(dic_V['63']-dic_V['72'])),
                        '3': (3./5.)*np.sqrt((39./14.)*np.pi)*(dic_V['42']-dic_V['53']) + (1./5.)*(np.sqrt((78./7.)*np.pi))*(-dic_V['71']),
                        '-3': (3./5.)*np.sqrt((39./14.)*np.pi)*(dic_V['43']+dic_V['52']) + (1./5.)*(np.sqrt((78./7.)*np.pi))*(dic_V['61']),
                        '4': (3./5.)*np.sqrt((26./7.)*np.pi)*(-dic_V['44']/2 + dic_V['55']/2) + np.sqrt((39./70.)*np.pi)*(-dic_V['62']+dic_V['73']),
                        '-4': (3./5.)*np.sqrt((26./7.)*np.pi)*(-dic_V['54']) + np.sqrt((39./70.)*np.pi)*(-dic_V['63']-dic_V['72']),
                        '5': (1./5.)*np.sqrt((429./14.)*np.pi)*(dic_V['64']-dic_V['75']),
                        '-5': (1./5.)*np.sqrt((429./14.)*np.pi)*(dic_V['65']+dic_V['74']),
                        '6':(1./5.)*(np.sqrt((429./7.)*np.pi))*(-dic_V['66']/2 + dic_V['77']/2),
                        '-6':(1./5.)*(np.sqrt((429./7.)*np.pi))*(-dic_V['76'])}
                    }

    else:
        print('ERROR: in from_Vint_to_Bkq')
        exit()

    #Adjust phases so the Bkq and Bkq' correspond to those of Goeller-Walrand
    for k in range(0,2*l+1,2):
        for q in range(0,k+1):
            dic_ckq[str(k)][str(q)] *= np.sqrt((2*k+1)/(4*np.pi))*(-1)**q
            if q!=0:
                dic_ckq[str(k)][str(-q)] *= -np.sqrt((2*k+1)/(4*np.pi))*(-1)**q

    return dic_ckq

def projection_LS(basis2, labels, basis1, bin=1e-4, J_label=False):   #!!!!!!!NOT USED (BUT WORKS)
    # calculates projection of basis2 on basis1 (the free ion case, with only e-e interaction)
    # expresses it in terms of SL basis or SLJ basis (if J_label = True)

    # print(labels)

    # LS_list = [labels[i][:2] for i in range(len(labels))]
    # LS = []
    # [LS.append(LS_list[i]) for i in range(len(LS_list)) if LS_list[i] not in LS]

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
            for key in labels[j+1].keys():   #dovrebbe essere solo 1
                stato = [key, matrix_coeff_re2[j,i]]
            states[i].append(stato)
            if J_label==True:
                key, value = stato[0], stato[1]
            else:
                key, value = stato[0][:2], stato[1]
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
        # print(tot)
        # print(matrix_coeff[:,i])   #combinazione lineare per basis2[:,i]
        # print(matrix_coeff_re2[:,i])
        # print(states[i])
        # print(states_red[i+1])
    # print(states_red)
    # print(states_red)
    # exit()
    return states_red
