import numpy as np
import scipy

#########CRYSTAL FIELD###############

def w_inp_dat(coord_tot, charges_tot, num=1, directory=None):
    """
    Creates a file named 'simpre.dat', which contains coordinates and charges, to be used as input for a SIMPRE calculation.
    As such, only the negative charges are considered.
    ----------
    Parameters:
    - coord_tot : 2darray
        cartesian coordinates of the ligands
    - charges_tot : sequence
        ligand charges (with their sign)
    - num : int
        SIMPRE index for a configuration
    """

    coord = []
    charges = []
    for i in range(coord_tot.shape[0]):
        if charges_tot[i]<=0:
            coord.append(coord_tot[i,:])
            charges.append(np.abs(charges_tot[i]))
    coord = np.array(coord)
    if coord.size==0:
        print('WARNING: Only positive charged ligands found. The "simpre.dat" file can not be generated.')
        return 0

    with open(directory+'simpre.dat', 'w') as f:
        f.write('1=Ce3+, 2=Pr3+, 3=Nd3+, 4=Pm3+, 5=Sm3+, 6=Tb3+, 7=Dy3+, 8=Ho3+, 9=Er3+, 10=Tm3+, 11=Yb3+, 12=user\n')
        f.write('\n')
        f.write(' '+f'{int(num):2.0f}'+'    !ion code (from 1 to 12, see above)\n')
        f.write(' '+f'{coord.shape[0]:2.0f}'+'    !number of effective charges\n')
        for i in range(len(coord[:,0])):
            f.write('{:-3}'.format(i+1))
            for j in range(3):
                f.write('{:13.7f}'.format(coord[i,j]))
            f.write('{:13.7f}'.format(charges[i])+'\n')
        f.write('\n')
        f.write('\n')
        f.close()

def read_data(filename, sph_flag = False):
    #the angles must be in radiants

    file = open(filename).readlines()
    coord_m = []
    charges = []
    labels = []
    rdisp = []
    for i,line in enumerate(file):
        splitline = line.split('\t')
        labels.append(splitline[0])
        charges.append(eval(splitline[1]))
        rdisp.append(float(splitline[2]))
        coord_m.append([float(splitline[j]) for j in range(3,len(splitline))])
    coord_m = np.array(coord_m)
    rdisp = np.array(rdisp)

    if not sph_flag:
        sph_coord = from_car_to_sph(coord_m)
    else:
        sph_coord = np.copy(coord_m)

    sph_coord[:,0] -= rdisp

    coord_m = from_sph_to_car(sph_coord)
    charges = np.array(charges)

    data = np.zeros((coord_m.shape[0],5), dtype='object')
    data[:,1:-1] = coord_m
    data[:,-1] = np.array(charges)
    data[:,0] = labels

    return data

def from_Vint_to_Bkq_2(l, dic_, reverse=False):
    #equivalent to from_Vint_to_Bkq but with the possibility of doing the reverse operation

    dic_conv = copy.deepcopy(dic_)

    if reverse:
        #fill the dictionary with zeros
        for k in range(0,2*l+1,2):
            if str(k) not in dic_conv.keys():
                dic_conv[str(k)] = {}
            for q in range(k,-k-1,-1):
                if str(q) not in dic_conv[str(k)].keys():
                    dic_conv[str(k)][str(q)] = 0.0
        V = []
        for k in range(0,2*l+1,2):
            for q in range(0,k+1):
                dic_conv[str(k)][str(q)] /= np.sqrt((2*k+1)/(4*np.pi))*(-1)**q
                V.append(dic_conv[str(k)][str(q)])
                if q!=0:
                    dic_conv[str(k)][str(-q)] /= np.sqrt((2*k+1)/(4*np.pi))*(-1)**(q+1)
                    V.append(dic_conv[str(k)][str(-q)])
    else:
        V = [dic_conv[key] for key in dic_conv.keys()]

    if l==2:

        B0_0 = [2./5.*np.sqrt(np.pi), 0, 2./5.*np.sqrt(np.pi), 0, 0, 2./5.*np.sqrt(np.pi), 0, 0, 0, 2./5.*np.sqrt(np.pi), 0, 0, 0, 0, 2./5.*np.sqrt(np.pi)]
        B2_0 = [np.sqrt(np.pi/5.)*(2.),0, np.sqrt(np.pi/5.),0, 0, np.sqrt(np.pi/5.), 0, 0, 0, np.sqrt(np.pi/5.)*(- 2.), 0, 0, 0, 0, np.sqrt(np.pi/5.)*(- 2.)]
        B2_1 = [0, 0, 0, -np.sqrt(4.*np.pi/5.)*1/np.sqrt(2), 0, 0, 0, -np.sqrt(4.*np.pi/5.)*( np.sqrt(3)/np.sqrt(2)), 0, 0, 0, 0, -np.sqrt(4.*np.pi/5.)*( np.sqrt(3)/np.sqrt(2)), 0, 0]
        B2_m1 = [0, np.sqrt(4.*np.pi/5.)*1/np.sqrt(2), 0, 0, 0, 0, 0, 0,  np.sqrt(4.*np.pi/5.)*(-np.sqrt(3)/np.sqrt(2)*(-1)), 0, 0, np.sqrt(4.*np.pi/5.)*(-np.sqrt(3)/np.sqrt(2)), 0, 0, 0]
        B2_2 = [0, 0, -np.sqrt(4.*np.pi/5.)*np.sqrt(3)/(2.*np.sqrt(2)), 0, 0, -np.sqrt(4.*np.pi/5.)*(np.sqrt(3)/(2.*np.sqrt(2))*(-1)), 0, 0, 0, 0, -np.sqrt(4.*np.pi/5.)*np.sqrt(2), 0, 0, 0, 0]
        B2_m2 = [0, 0, 0, 0, -np.sqrt(4.*np.pi/5.)* np.sqrt(3)/np.sqrt(2), 0, -np.sqrt(4.*np.pi/5.)*(-np.sqrt(2)), 0, 0, 0, 0, 0, 0, 0, 0]
        B4_0 = [np.sqrt(np.pi)/5.*(6.),0, np.sqrt(np.pi)/5.*(-4.),0, 0, np.sqrt(np.pi)/5.*(-4.), 0, 0, 0, np.sqrt(np.pi)/5., 0, 0, 0, 0, np.sqrt(np.pi)/5.]
        B4_1 = [0, 0, 0, 2*np.sqrt(2.*np.pi/5.)*(-np.sqrt(3)/np.sqrt(2)), 0, 0, 0, 2*np.sqrt(2.*np.pi/5.)*1./(2.*np.sqrt(2)), 0, 0, 0, 0, 2*np.sqrt(2.*np.pi/5.)*(1./(2.*np.sqrt(2))), 0, 0]
        B4_m1 = [0, 2*np.sqrt(2.*np.pi/5.)*(np.sqrt(3)/np.sqrt(2)), 0, 0, 0, 0, 0, 0, 2*np.sqrt(2.*np.pi/5.)*(1./(2.*np.sqrt(2))*(-1)), 0, 0, 2*np.sqrt(2.*np.pi/5.)*1./(2.*np.sqrt(2)), 0, 0, 0]
        B4_2 = [0, 0, 2.*np.sqrt(2.*np.pi/5.)*(-1)/2., 0, 0, 2.*np.sqrt(2.*np.pi/5.)*1/2, 0, 0, 0, 0, 2.*np.sqrt(2.*np.pi/5.)*(np.sqrt(3)/2.), 0, 0, 0, 0]
        B4_m2 = [0, 0, 0, 0, 2.*np.sqrt(2.*np.pi/5.)*(-1), 0, 2.*np.sqrt(2.*np.pi/5.)*(-np.sqrt(3)/2.), 0, 0, 0, 0, 0, 0, 0, 0]
        B4_3 = [0, 0, 0, 0, 0, 0, 0, np.sqrt(7*np.pi/5.), 0, 0, 0, 0, np.sqrt(7*np.pi/5.)*(-1), 0, 0]
        B4_m3 = [0, 0, 0, 0, 0, 0, 0, 0, np.sqrt(7*np.pi/5.), 0, 0, np.sqrt(7*np.pi/5.), 0, 0, 0]
        B4_4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, np.sqrt(7.*np.pi/10.)*(-1), 0, 0, 0, 0, np.sqrt(7.*np.pi/10.)]
        B4_m4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -np.sqrt(7.*np.pi/5.)*np.sqrt(2), 0]

        conv_matrix = np.vstack((B0_0, B2_0, B2_1, B2_m1, B2_2, B2_m2, B4_0, B4_1, B4_m1, B4_2, B4_m2, B4_3, B4_m3, B4_4, B4_m4))


    elif l==3:

        B0_0 = [(2./7.)*np.sqrt(np.pi), 0, (2./7.)*np.sqrt(np.pi), 0, 0, (2./7.)*np.sqrt(np.pi), 0, 0, 0, (2./7.)*np.sqrt(np.pi), 0, 0, 0, 0, (2./7.)*np.sqrt(np.pi), 0, 0, 0, 0, 0, (2./7.)*np.sqrt(np.pi), 0, 0, 0, 0, 0, 0, (2./7.)*np.sqrt(np.pi)]
        B2_0 = [(2./7.)*np.sqrt(5*np.pi), 0, (3./14.)*np.sqrt(5*np.pi), 0, 0, (3./14.)*np.sqrt(5*np.pi), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, - (5./14.)*np.sqrt(5*np.pi), 0, 0, 0, 0, 0, 0, - (5./14.)*np.sqrt(5*np.pi)]
        B2_1 = [0, 0, 0, - (1./7.)*np.sqrt(5*np.pi), 0, 0, 0,  (5./14.)*np.sqrt(3*np.pi)*(-1), 0, 0, 0, 0, (5./14.)*np.sqrt(3*np.pi)*(-1), 0, 0, 0, 0, 0, (5./14.)*np.sqrt(5*np.pi)*(-1), 0, 0, 0, 0, 0, 0, (5./14.)*np.sqrt(5*np.pi)*(-1), 0, 0]
        B2_m1 = [0, (1./7.)*np.sqrt(5*np.pi), 0, 0, 0, 0, 0, 0, (5./14.)*np.sqrt(3*np.pi), 0, 0, (5./14.)*np.sqrt(3*np.pi)*(-1), 0, 0, 0, 0, 0, 0, 0, (5./14.)*np.sqrt(5*np.pi), 0, 0, 0, 0, (5./14.)*np.sqrt(5*np.pi)*(-1), 0, 0, 0]
        B2_2 = [0, 0, (1./7.)*np.sqrt(30*np.pi)*(-1/2), 0, 0, (1./7.)*np.sqrt(30*np.pi)*(1/2), 0, 0, 0, 0, (5./7.)*np.sqrt(2*np.pi)*(-1), 0, 0, 0, 0, 0, (5./7.)*np.sqrt(np.pi/2)*(-1), 0, 0, 0, 0, 0, 0, (5./7.)*np.sqrt(np.pi/2)*(-1), 0, 0, 0, 0]
        B2_m2 = [0, 0, 0, 0, (1./7.)*np.sqrt(30*np.pi)*(-1), 0, (5./7.)*np.sqrt(2*np.pi), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (5./7.)*np.sqrt(np.pi/2), 0, 0, 0, 0, (5./7.)*np.sqrt(np.pi/2)*(-1), 0, 0, 0, 0, 0]
        B4_0 = [(1./7.)*np.sqrt(np.pi)*(6), 0, (1./7.)*np.sqrt(np.pi), 0, 0, (1./7.)*np.sqrt(np.pi), 0, 0, 0, -np.sqrt(np.pi), 0, 0, 0, 0,  -np.sqrt(np.pi), 0, 0, 0, 0, 0, (1./7.)*np.sqrt(np.pi)*(3), 0, 0, 0, 0, 0, 0, (1./7.)*np.sqrt(np.pi)*(3)]
        B4_1 = [0, 0, 0, (1./7.)*np.sqrt(30*np.pi)*(-1), 0, 0, 0, (4./7.)*np.sqrt(2*np.pi)*(-1), 0, 0, 0, 0, (4./7.)*np.sqrt(2*np.pi)*(-1), 0, 0, 0, 0, 0, (1./7.)*np.sqrt(30*np.pi), 0, 0, 0, 0, 0, 0, (1./7.)*np.sqrt(30*np.pi), 0, 0]
        B4_m1 = [0, (1./7.)*np.sqrt(30*np.pi), 0, 0, 0, 0, 0, 0, (4./7.)*np.sqrt(2*np.pi), 0, 0, (4./7.)*np.sqrt(2*np.pi)*(-1), 0, 0, 0, 0, 0, 0, 0, (1./7.)*np.sqrt(30*np.pi)*(-1), 0, 0, 0, 0, (1./7.)*np.sqrt(30*np.pi), 0, 0, 0]
        B4_2 = [0, 0, (2./7.)*np.sqrt(10*np.pi)*(-1/2), 0, 0, (2./7.)*np.sqrt(10*np.pi)*(1/2), 0, 0, 0, 0, (1./7.)*np.sqrt(6*np.pi)*(-1), 0, 0, 0, 0, 0, (3./7.)*np.sqrt(6*np.pi), 0, 0, 0, 0, 0, 0, (3./7.)*np.sqrt(6*np.pi), 0, 0, 0, 0]
        B4_m2 = [0, 0, 0, 0, (2./7.)*np.sqrt(10*np.pi)*(-1), 0, (1./7.)*np.sqrt(6*np.pi), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (3./7.)*np.sqrt(6*np.pi)*(-1), 0, 0, 0, 0, (3./7.)*np.sqrt(6*np.pi), 0, 0, 0, 0, 0]
        B4_3 = [0, 0, 0, 0, 0, 0, 0, np.sqrt((2./7.)*np.pi), 0, 0, 0, 0, np.sqrt((2./7.)*np.pi)*(-1), 0, 0, 0, 0, 0, 0, 0, 0, 3*np.sqrt((2./7.)*np.pi), 0, 0, 0, 0, 0, 0]
        B4_m3 = [0, 0, 0, 0, 0, 0, 0, 0, np.sqrt((2./7.)*np.pi), 0, 0, np.sqrt((2./7.)*np.pi), 0, 0, 0, - 3*np.sqrt((2./7.)*np.pi), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        B4_4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, np.sqrt((10./7.)*np.pi)*(-1/2), 0, 0, 0, 0, np.sqrt((10./7.)*np.pi)*(1/2), 0, np.sqrt((6./7.)*np.pi), 0, 0, 0, 0, 0, 0, np.sqrt((6./7.)*np.pi)*(-1), 0, 0, 0, 0] 
        B4_m4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, np.sqrt((10./7.)*np.pi)*(-1), 0, 0, 0, np.sqrt((6./7.)*np.pi), 0, 0, 0, 0, np.sqrt((6./7.)*np.pi), 0, 0, 0, 0, 0]
        B6_0 = [(1./7.)*np.sqrt(13*np.pi)*(2), 0, (1./7.)*np.sqrt(13*np.pi)*(-(3./2.)), 0, 0, (1./7.)*np.sqrt(13*np.pi)*(-(3./2.)), 0, 0, 0, (1./7.)*np.sqrt(13*np.pi)*(3./5.), 0, 0, 0, 0, (1./7.)*np.sqrt(13*np.pi)*(3./5.), 0, 0, 0, 0, 0, -(1./7.)*np.sqrt(13*np.pi)*(1./10.), 0, 0, 0, 0, 0, 0, -(1./7.)*np.sqrt(13*np.pi)*(1./10.)]
        B6_1 = [0, 0, 0, np.sqrt((13./7.)*np.pi)*(-1), 0, 0, 0, np.sqrt((13./7.)*np.pi)*(1./2.)*np.sqrt(3./5.), 0, 0, 0, 0, np.sqrt((13./7.)*np.pi)*(1./2.)*np.sqrt(3./5.), 0, 0, 0, 0, 0, np.sqrt((13./7.)*np.pi)*(1./10.)*(-1), 0, 0, 0, 0, 0, 0, np.sqrt((13./7.)*np.pi)*(1./10.)*(-1), 0, 0]
        B6_m1 = [0, np.sqrt((13./7.)*np.pi), 0, 0, 0, 0, 0, 0, np.sqrt((13./7.)*np.pi)*(1./2.)*np.sqrt(3./5.)*(-1), 0, 0, np.sqrt((13./7.)*np.pi)*(1./2.)*np.sqrt(3./5.), 0, 0, 0, 0, 0, 0, 0, np.sqrt((13./7.)*np.pi)*(1./10.), 0, 0, 0, 0, np.sqrt((13./7.)*np.pi)*(1./10.)*(-1), 0, 0, 0]
        B6_2 = [0, 0, np.sqrt((13./7.)*np.pi)*(-1./2.)*np.sqrt(3./5.), 0, 0, np.sqrt((13./7.)*np.pi)*(1./2.)*np.sqrt(3./5.), 0, 0, 0, 0, np.sqrt((13./7.)*np.pi)*4/5, 0, 0, 0, 0, 0, np.sqrt((13./7.)*np.pi)*(1./5.)*(-1), 0, 0, 0, 0, 0, 0, np.sqrt((13./7.)*np.pi)*(1./5.)*(-1), 0, 0, 0, 0]
        B6_m2 = [0, 0, 0, 0, np.sqrt((13./7.)*np.pi)*(-np.sqrt(3./5.)), 0, np.sqrt((13./7.)*np.pi)*(- 4/5), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, np.sqrt((13./7.)*np.pi)*(1./5.), 0, 0, 0, 0, np.sqrt((13./7.)*np.pi)*(1./5.)*(-1), 0, 0, 0, 0, 0]
        B6_3 = [0, 0, 0, 0, 0, 0, 0, (3./5.)*np.sqrt((39./14.)*np.pi), 0, 0, 0, 0, (3./5.)*np.sqrt((39./14.)*np.pi)*(-1), 0, 0, 0, 0, 0, 0, 0, 0, (1./5.)*(np.sqrt((78./7.)*np.pi))*(-1), 0, 0, 0, 0, 0, 0]
        B6_m3 = [0, 0, 0, 0, 0, 0, 0, 0, (3./5.)*np.sqrt((39./14.)*np.pi), 0, 0, (3./5.)*np.sqrt((39./14.)*np.pi), 0, 0, 0, (1./5.)*(np.sqrt((78./7.)*np.pi)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        B6_4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, (3./5.)*np.sqrt((26./7.)*np.pi)*(-1/2), 0, 0, 0, 0, (3./5.)*np.sqrt((26./7.)*np.pi)*(1/2), 0, np.sqrt((39./70.)*np.pi)*(-1), 0, 0, 0, 0, 0, 0, np.sqrt((39./70.)*np.pi), 0, 0, 0, 0]
        B6_m4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (3./5.)*np.sqrt((26./7.)*np.pi)*(-1), 0, 0, 0, np.sqrt((39./70.)*np.pi)*(-1), 0, 0, 0, 0, np.sqrt((39./70.)*np.pi)*(-1), 0, 0, 0, 0, 0]
        B6_5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1./5.)*np.sqrt((429./14.)*np.pi), 0, 0, 0, 0, 0, 0, (1./5.)*np.sqrt((429./14.)*np.pi)*(-1), 0, 0]
        B6_m5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1./5.)*np.sqrt((429./14.)*np.pi), 0, 0, 0, 0, (1./5.)*np.sqrt((429./14.)*np.pi), 0, 0, 0]
        B6_6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1./5.)*(np.sqrt((429./7.)*np.pi))*(-1/2), 0, 0, 0, 0, 0, 0, (1./5.)*(np.sqrt((429./7.)*np.pi))*(1/2)]  
        B6_m6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (1./5.)*(np.sqrt((429./7.)*np.pi))*(-1), 0] 

        conv_matrix = np.vstack((B0_0, B2_0, B2_1, B2_m1, B2_2, B2_m2, B4_0, B4_1, B4_m1, B4_2, B4_m2, B4_3, B4_m3, B4_4, B4_m4, B6_0, B6_1, B6_m1, B6_2, B6_m2, B6_3, B6_m3, B6_4, B6_m4, B6_5, B6_m5, B6_6, B6_m6))

    
    if not reverse:
        B = np.dot(conv_matrix, V)
        dic_Bkq = {}
        count = 0
        for k in range(0,2*l+1,2):
            dic_Bkq[str(k)] = {}
            for q in range(0,k+1):
                dic_Bkq[str(k)][str(q)] = np.sqrt((2*k+1)/(4*np.pi))*(-1)**q*B[count]    # to match the "reverse" phase convention ([+q] + (-1)**q [-q])  ...
                count += 1
                if q!=0:
                    dic_Bkq[str(k)][str(-q)] = np.sqrt((2*k+1)/(4*np.pi))*(-1)**(q+1)*B[count]    # to match the "reverse" phase convention ([+q] + (-1)**q [-q])  ...
                    count += 1 

        return dic_Bkq
    
    else:
        B = np.dot(np.linalg.inv(conv_matrix), V)
        dic_V = {}
        count = 0
        for k in range(0,2*l+1):
            for q in range(0,k+1):
                dic_V[str(k+1)+str(q+1)] = B[count]
                count += 1

        return dic_V

def from_AOM_to_Vint(dic_AOM, conf):
    #conversion from AOM dict [es, eps, epc, theta, phi, chi] to <|V|> = sum(lig) sum(modes) D*D*e

    if conf[0]=='d':
        l=2
    else:
        l=3

    ligand_label = [key for key in dic_AOM.keys()]
    matrix_elements = [value for value in dic_AOM.values()]
    matrix_elements = np.array(matrix_elements)
    ee = matrix_elements[:,:3]    # e_sigma, e_pi_s, e_pi_c
    theta = matrix_elements[:,3]
    phi = matrix_elements[:,4]
    if matrix_elements.shape[1]==5:
        chi = np.zeros(matrix_elements.shape[0])
    else:
        chi = matrix_elements[:,-1]

    D = np.zeros((2*l+1,3))
    V = np.zeros((2*l+1,2*l+1))
    for i in range(matrix_elements.shape[0]): #per ogni legante
        theta[i] *= np.pi/180
        phi[i] *= np.pi/180
        chi[i] *= np.pi/180

        a1 = np.cos(phi[i])*np.cos(theta[i])*np.cos(chi[i])-np.sin(phi[i])*np.sin(chi[i])
        a2 = np.sin(phi[i])*np.cos(theta[i])*np.cos(chi[i])+np.cos(phi[i])*np.sin(chi[i])
        a3 = -np.sin(theta[i])*np.cos(chi[i])
        b1 = -np.cos(phi[i])*np.cos(theta[i])*np.sin(chi[i])-np.sin(phi[i])*np.cos(chi[i])
        b2 = -np.sin(phi[i])*np.cos(theta[i])*np.sin(chi[i])+np.cos(phi[i])*np.cos(chi[i])
        b3 = np.sin(theta[i])*np.sin(chi[i])
        g1 = np.cos(phi[i])*np.sin(theta[i])
        g2 = np.sin(phi[i])*np.sin(theta[i])
        g3 = np.cos(theta[i])

        if l==3:
            #Tab 1 Urland 1976 Chem. Phys. 14, 393-401
            D[0,0] = 0.5*g3*(5.*g3**2-3)
            D[0,1] = np.sqrt(3./2.)*b3*0.5*(5*g3**2-1)
            D[0,2] = np.sqrt(3./2.)*a3*0.5*(5*g3**2-1)
            D[1,0] = (1./4.)*np.sqrt(6)*g2*(5*g3**2-1)
            D[1,1] = (1./4.)*b2*(5*g3**2-1)+(5./2.)*g2*g3*b3
            D[1,2] = (1./4.)*a2*(5*g3**2-1)+(5./2.)*g2*g3*a3
            D[2,0] = (1./4.)*np.sqrt(6)*g1*(5*g3**2-1)
            D[2,1] = (1./4.)*b1*(5*g3**2-1)+(5./2.)*g1*g3*b3
            D[2,2] = (1./4.)*a1*(5*g3**2-1)+(5./2.)*g1*g3*a3
            D[3,0] = np.sqrt(15)*g1*g2*g3
            D[3,1] = np.sqrt(5./2.)*(b1*g2*g3+b2*g1*g3+b3*g1*g2)
            D[3,2] = np.sqrt(5./2.)*(a1*g2*g3+a2*g1*g3+a3*g1*g2)
            D[4,0] = 0.5*np.sqrt(15)*g3*(g1**2-g2**2)
            D[4,1] = np.sqrt(5./8.)*(g1**2*b3-g2**2*b3+2*b1*g1*g3-2*b2*g2*g3)
            D[4,2] = np.sqrt(5./8.)*(g1**2*a3-g2**2*a3+2*a1*g1*g3-2*a2*g2*g3)
            D[5,0] = (1./4.)*np.sqrt(10)*g2*(3*g1**2-g2**2)
            D[5,1] = (1./4.)*np.sqrt(15)*b2*(g1**2-g2**2)+0.5*np.sqrt(15)*g1*g2*b1
            D[5,2] = (1./4.)*np.sqrt(15)*a2*(g1**2-g2**2)+0.5*np.sqrt(15)*g1*g2*a1
            D[6,0] = (1./4.)*np.sqrt(10)*g1*(g1**2-3*g2**2)
            D[6,1] = (1./4.)*np.sqrt(15)*b1*(g1**2-g2**2)-0.5*np.sqrt(15)*g1*g2*b2
            D[6,2] = (1./4.)*np.sqrt(15)*a1*(g1**2-g2**2)-0.5*np.sqrt(15)*g1*g2*a2
        elif l==2:
            #eq 9 Gerloch & McMeeking 1975
            D[0,0] = 0.5*(3*g3**2-1)
            D[0,1] = np.sqrt(3)*a3*g3
            D[0,2] = np.sqrt(3)*b3*g3
            D[1,0] = np.sqrt(3)*g1*g3
            D[1,1] = a1*g3+a3*g1
            D[1,2] = b1*g3+b3*g1
            D[2,0] = np.sqrt(3)*g2*g3
            D[2,1] = a2*g3+a3*g2
            D[2,2] = b2*g3+b3*g2
            D[3,0] = np.sqrt(3)*g1*g2
            D[3,1] = a1*g2+a2*g1
            D[3,2] = b1*g2+b2*g1
            D[4,0] = np.sqrt(3)*0.5*(g1**2-g2**2)
            D[4,1] = a1*g1-a2*g2
            D[4,2] = b1*g1-b2*g2
        else:
            print('ERROR: l != 2 and l != 3 in from_AOM_to_Vint')

        #M = A*e  (eq 10 Gerloch & McMeeking 1975)
        for ii in range(0,2*l+1):
            for jj in range(0,2*l+1):
                V[ii,jj] += ee[i,0]*D[ii,0]*D[jj,0] + ee[i,1]*D[ii,2]*D[jj,2] + ee[i,2]*D[ii,1]*D[jj,1]  #in input li fornisco come x-y ma D è costruita per y-x

    dic_V = {}
    for i in range(1,2*l+1+1):
        for j in range(1,i+1):
            dic_V['{}{}'.format(i,j)] = V[i-1,j-1]

    return dic_V

def from_Aqkrk_to_Bkq(Aqkrk, revers=False):  
    # conversion from Akq<rk> of stevens to Bkq of wybourne (o Lkq according to: https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/manual/node130.html)

    dic_bkq = {}
    for key1 in Aqkrk.keys():
        dic_bkq[key1] = {}
        for key2, value in Aqkrk[key1].items():
            if revers==False:
                dic_bkq[key1].update({key2: value*conv_Aqkrk_bkq(int(key1),np.abs(int(key2)))})
            else:
                dic_bkq[key1].update({key2: value/conv_Aqkrk_bkq(int(key1),np.abs(int(key2)))})

    return dic_bkq

def sph_harm(l,m,theta,phi):
    yL = scipy.special.lpmn(m, l, np.cos(theta))[0][-1][-1]
    y = np.sqrt(fact(l-m)/fact(l+m))*yL*np.exp(1j*m*phi)
    return y

def calc_Bqk(data, conf, sph_flag = False, sth_param = False, bin=1e-9):
    # calc Stevens coefficients, B^q_k, from data in the point charge model

    if conf[0]=='d':
        l=2
    else:
        l=3

    dic_Aqk_rk = calc_Aqkrk(data, conf, sph_flag, sth_param, bin)
    dic_Bqk = {}
    for k in range(2,2*l+1,2):
        dic_Bqk[f'{k}'] = {}
        for q in range(-k,k+1):
            dic_Bqk[f'{k}'].update({f'{q}': dic_Aqk_rk[f'{k}'][f'{q}']*Stev_coeff(str(k), conf)})

    return dic_Bqk

def calc_Aqkrk(data, conf, sph_flag = False, sth_param = False, bin=1e-9):
    #in data is (N. ligands x 5)
    #[labels, x, y, z, charge] if sph_falg==False
    #the calculation is performed in the hard point charge model
    #equation 9a, 9b, 9c from Software package SIMPRE - Revisited (M. Karbowiak and C. Rudowicz)

    import scipy.special

    au_conv = [scipy.constants.physical_constants['hartree-inverse meter relationship'][0]*1e-2, 1.889725989] #convertion from atomic unit

    if conf[0]=='d':
        l=2
    else:
        l=3

    if sph_flag==False:
        coord_car = data[:,1:-1]
        coord_sph = from_car_to_sph(coord_car)
    else:
        coord_sph = data[:,1:-1]
        coord_sph[:,1:] *= np.pi/180

    Aqkrk = {}
    for k in range(2,2*l+1,2):
        Aqkrk[f'{k}'] = {}
        r_val = r_expect(str(k), conf)
        for q in range(-k,k+1):
            pref = plm(k,np.abs(q))*(4*np.pi/(2*k+1))
            somma = 0
            for i in range(data.shape[0]):
                sphharmp = scipy.special.sph_harm(np.abs(q), k, coord_sph[i,2],coord_sph[i,1])  #sph_harm(k,q,coord_sph[i,1], coord_sph[i,2])/np.sqrt(4*np.pi/(2*k+1))
                sphharmm = scipy.special.sph_harm(-np.abs(q), k, coord_sph[i,2],coord_sph[i,1])  #sph_harm(k,-q,coord_sph[i,1], coord_sph[i,2])/np.sqrt(4*np.pi/(2*k+1))
                r = coord_sph[i,0]*au_conv[1]
                if q==0:
                    somma += pref*(data[i,-1]*au_conv[0]/r**(k+1))*scipy.special.sph_harm(0, k, coord_sph[i,2],coord_sph[i,1]).real
                elif q>0:
                    somma += pref*(data[i,-1]*au_conv[0]/r**(k+1))*(1/np.sqrt(2))*(sphharmm + (-1)**q*sphharmp).real
                elif q<0:
                    somma += (1j*pref*(data[i,-1]*au_conv[0]/r**(k+1))*(1/np.sqrt(2))*(sphharmm - (-1)**q*sphharmp)).real

            if sth_param == True:
                value = (1-sigma_k(str(k), conf))*somma*r_val
            else:
                value = somma*r_val
            if np.abs(value)<=bin:
                Aqkrk[f'{k}'].update({f'{q}': 0})
            else:
                Aqkrk[f'{k}'].update({f'{q}': value})

    return Aqkrk

def calc_Bkq(data, conf, sph_flag = False, sth_param = False, bin=1e-9):  
    #eq 11a-11b-11c-12 from Software package SIMPRE - Revisited (M. Karbowiak and C. Rudowicz)

    import scipy.special

    au_conv = [scipy.constants.physical_constants['hartree-inverse meter relationship'][0]*1e-2, 1.889725989] #convertion from atomic unit

    if conf[0]=='d':
        l=2
    else:
        l=3

    if sph_flag==False:
        coord_car = data[:,1:-1]
        coord_sph = from_car_to_sph(coord_car)
    else:
        coord_sph = data[:,1:-1]
        coord_sph[:,1:] *= np.pi/180

    Bkq = {}
    for k in range(2,2*l+1,2):
        Bkq[f'{k}'] = {}
        prefac = np.sqrt(4*np.pi/(2*k+1))
        r_val = r_expect(str(k), conf)
        for q in range(0,k+1):
            somma = 0
            for i in range(data.shape[0]):
                r = coord_sph[i,0]*au_conv[1]
                sphharm = scipy.special.sph_harm(q, k, coord_sph[i,2],coord_sph[i,1])
                if q==0:
                    somma += prefac*(data[i,-1]*au_conv[0]/r**(k+1))*sphharm.real
                else:
                    somma += (-1)**q*prefac*(data[i,-1]*au_conv[0]/r**(k+1))*sphharm
            if sth_param == True:
                value = (1-sigma_k(str(k), conf))*somma*r_val
            else:
                value = somma*r_val
            if np.abs(value)<=bin:
                Bkq[f'{k}'].update({f'{q}': 0})
            else:
                if q!=0:
                    Bkq[f'{k}'].update({f'{q}': value.real})
                    Bkq[f'{k}'].update({f'{-q}': value.imag})
                else:
                    Bkq[f'{k}'].update({f'{q}': value})

    return Bkq


def rota_LF(l, dic_Bkq, A=0, B=0, C=0):
    #rotation of Bkq with Euler angles
    #the angles must be in radiants
    #convention of D-matrix: Z-Y-Z

    dic_Bkq_new = {}

    for k in range(2,2*l+1,2):
        D = np.zeros((k*2+1,k*2+1), dtype='complex128')
        for ii,m1 in enumerate(range(-k,k+1)):
            for jj,m in enumerate(range(-k,k+1)):
                D[ii,jj] = Wigner_coeff.Wigner_Dmatrix(k, m1, m, A, B, C)

        dic_Bkq_new[str(k)] = {}

        Bkq_vec = []
        for q in range(-k,k+1):
            if q!=0:
                if q<0:
                    Bkq_vec.append(dic_Bkq[str(k)][str(np.abs(q))]+1j*dic_Bkq[str(k)][str(-np.abs(q))])
                else:
                    if q%2!=0:
                        Bkq_vec.append(-dic_Bkq[str(k)][str(np.abs(q))]+1j*dic_Bkq[str(k)][str(-np.abs(q))])
                    else:
                        Bkq_vec.append(dic_Bkq[str(k)][str(np.abs(q))]-1j*dic_Bkq[str(k)][str(-np.abs(q))])
            else:
                Bkq_vec.append(dic_Bkq[str(k)][str(np.abs(q))]+1j*0)


        Bkq_vec = np.array(Bkq_vec)
        Bkq_vec_new = D@Bkq_vec
        dic_Bkq_new[str(k)][str(0)] = Bkq_vec_new[k].real

        for i,q in enumerate(range(k,0,-1)):
            dic_Bkq_new[str(k)][str(q)] = Bkq_vec_new[i].real
            dic_Bkq_new[str(k)][str(-q)] = Bkq_vec_new[i].imag

    return dic_Bkq_new

#@cron
def rota_LF_quat(l, dic_Bkq, R, dict=None, coeff=None):
    #rotation of Bkq with quaternions

    dic_Bkq_new = {}

    for k in range(2,2*l+1,2):

        D = Wigner_coeff.Wigner_Dmatrix_quat_complete(k, R, dict = dict, coeff=coeff)

        dic_Bkq_new[str(k)] = {}

        Bkq_vec = []
        for q in range(-k,k+1):
            if q!=0:
                if q<0:
                    Bkq_vec.append(dic_Bkq[str(k)][str(np.abs(q))]+1j*dic_Bkq[str(k)][str(-np.abs(q))])
                else:
                    if q%2!=0:
                        Bkq_vec.append(-dic_Bkq[str(k)][str(np.abs(q))]+1j*dic_Bkq[str(k)][str(-np.abs(q))])
                    else:
                        Bkq_vec.append(dic_Bkq[str(k)][str(np.abs(q))]-1j*dic_Bkq[str(k)][str(-np.abs(q))])
            else:
                Bkq_vec.append(dic_Bkq[str(k)][str(np.abs(q))]+1j*0)


        Bkq_vec = np.array(Bkq_vec)
        Bkq_vec_new = D@Bkq_vec
        dic_Bkq_new[str(k)][str(0)] = Bkq_vec_new[k].real

        for i,q in enumerate(range(k,0,-1)):
            dic_Bkq_new[str(k)][str(q)] = Bkq_vec_new[i].real
            dic_Bkq_new[str(k)][str(-q)] = Bkq_vec_new[i].imag

    return dic_Bkq_new

def read_DWigner_quat():

    def coeff_op(k):
        matrix = np.ones((2*k+1,2*k+1))
        for i,ii in enumerate(range(k,-k-1,-1)):
            for j,jj in enumerate(range(k,-k-1,-1)):
                if j==0 and i==0:
                    pass
                elif j>0:
                    matrix[i,j] = matrix[i,j-1]*np.sqrt(k*(k+1)-(jj+1)*((jj+1)-1))
                    #print(j,jj,np.sqrt(k*(k+1)-(jj+1)*((jj+1)-1)))
            if i>0:
                matrix[i,:] = matrix[i-1,:]*np.sqrt(k*(k+1)-(ii+1)*((ii+1)-1))
        return matrix

    filename = ['tables/tab_wignerDquat.txt', 'tables/tab_wignerDquat_coeff_t.txt']
    list_dict = []
    for ii in range(len(filename)):
        file = open(filename[ii]).readlines()

        dict = {}
        last = 0
        for i,line in enumerate(file):
            if 'L=' in line:
                line = line.replace('\n','')
                splitline = line.split('=')
                dict[int(splitline[-1])]={}
                last = int(splitline[-1])
            elif 'L=' not in line and line!='\n':
                line = line.replace('\n','')
                splitline=line.split(' ')
                dict[last].update({splitline[0]+':'+splitline[1]:splitline[-1]})
        list_dict.append(dict)

    dic_matrix_coeff = {}
    for k in list_dict[0].keys():
        matrix_coeff = np.zeros((2*k+1,2*k+1), dtype='object')
        matrix_divide = coeff_op(k)
        for i,ii in enumerate(range(k,-1,-1)):
            #print(list(range(ii,-1,-1)))
            for j,jj in enumerate(range(ii,-1,-1)):
                j += i
                key = str(ii)+':'+str(jj)
                #print(key)
                #print(key, i,j, ii, jj, list_dict[1][k][key], matrix_divide[i,j])
                matrix_coeff[i,j] = np.abs(eval(list_dict[1][k][key]))/matrix_divide[i,j] #str(ii)+':'+str(jj)
                matrix_coeff[i,-j-1] = matrix_coeff[i,j]
                matrix_coeff[-i-1,j] = matrix_coeff[i,j]
                matrix_coeff[-i-1,-j-1] = matrix_coeff[i,j]
                matrix_coeff[j,i] = matrix_coeff[i,j]
                matrix_coeff[j,-i-1] = matrix_coeff[i,j]
                matrix_coeff[-j-1,i] = matrix_coeff[i,j]
                matrix_coeff[-j-1,-i-1] = matrix_coeff[i,j]
        #print(matrix_coeff)
        # # print(matrix_divide)
        #exit()
        dic_matrix_coeff[k]=matrix_coeff
    # dic_matrix_coeff = {}
    # for k in list_dict[0].keys():
    #     matrix_coeff = np.zeros((2*k+1,2*k+1))
    #     for key,value in list_dict[1][k].items():
    #         idx = [int(ii) for ii in key.split(':')]
    #         i = np.abs(idx[0]-k)
    #         j = np.abs(idx[1]-k)
    #         matrix_coeff[i,j] = np.abs(eval(value))
    #     #matrix_coeff /= coeff_op(k)
    #
    #     print(coeff_op(k))
    #     exit()
    #     dic_matrix_coeff[k]=matrix_coeff/coeff_op(k)
    # pprint(list_dict[0])
    # exit()

    return list_dict[0], dic_matrix_coeff

def R_zyz(alpha, beta, gamma):
    #rotation matrix zyz (active rotation)

    A = alpha*np.pi/180
    B = beta*np.pi/180
    C = gamma*np.pi/180

    ca = np.cos(A)
    sa = np.sin(A)
    cb = np.cos(B)
    sb = np.sin(B)
    cg = np.cos(C)
    sg = np.sin(C)

    R = np.array([[ca*cb*cg-sa*sg,-cg*sa-ca*cb*sg,ca*sb],[ca*sg+cb*cg*sa,ca*cg-cb*sa*sg,sa*sb],[-cg*sb,sb*sg,cb]])

    return R

def from_matrix_to_result(matrix):
    w,v = diagonalisation(matrix)  
    result = np.vstack((w,v))
    result = np.copy(result[:, result[0,:].argsort()])
    return result

def order_two_level_dict(dict1):
    dict_or = {}
    key1 = np.sort(np.array([eval(i) for i in dict1.keys()]))
    key2 = []
    for key in key1:
        key2.append(np.sort(np.array([eval(i) for i in dict1[str(key)].keys()])))
    for i,keyi in enumerate(key1):
        dict_or[str(keyi)] = {}
        for j,keyj in enumerate(key2[i]):
            try:
                dict_or[str(keyi)].update({str(keyj):dict1[str(keyi)][str(keyj)]})
            except:
                dict_or[str(keyi)] = {}
                dict_or[str(keyi)].update({str(keyj):dict1[str(keyi)][str(keyj)]})
    return dict_or

def Freeion_charge_dist(theta, phi, A2, A4, A6, r=1, bin=1e-10):
    #Graphical representation and Discussion of the Charge Density
    #from Ch. III p 292 of Sievers "Asphericity of 4f-Shells in their Hund's rule ground states" (1981)
    #in scipy.special.sph_harm the angles are defined in the opposite way

    c2 = A2/np.sqrt(4*np.pi/(2*2+1))
    c4 = A4/np.sqrt(4*np.pi/(2*4+1))
    c6 = A6/np.sqrt(4*np.pi/(2*6+1))

    val = 3/(4*np.pi) + c2*scipy.special.sph_harm(0,2,phi,theta).real + c4*scipy.special.sph_harm(0,4,phi,theta).real + c6*scipy.special.sph_harm(0,6,phi,theta).real
    if np.abs(val)<bin:
        val=0
    return (val)**(1/3)

def coeff_multipole_moments(conf, J, M, L=0, S=0):
    #based on calculation of charge density based on the Wigner-Eckart Theorem
    #from Ch. II p 290 of Sievers "Asphericity of 4f-Shells in their Hund's rule ground states" (1981)
    #Ak = sqrt(4*pi/(2*k+1))ck

    coeff = []
    for k in range(2,6+1,2):
        if int(conf[1:])>7:
            pref = (-1)**(J-M)*7/np.sqrt(4*np.pi)*Wigner_coeff.threej_symbol([[J,k,J],[-M,0,M]])/Wigner_coeff.threej_symbol([[J,k,J],[-J,0,J]])
            pref2 = np.sqrt(2*k+1)*Wigner_coeff.threej_symbol([[k,3,3],[0,0,0]])
            if k==0:  #not active, just for completness
                delta = 1
            else:
                delta = 0
            somma = 0
            for i in range(1,int(conf[1:])-7+1):
                somma += (-1)**i*Wigner_coeff.threej_symbol([[k,3,3],[0,4-i,i-4]])
            value = pref*(somma*pref2+delta)*np.sqrt(4*np.pi/(2*k+1))
            coeff.append(value)
        else:
            pref = (-1)**(2*J-M+L+S)*7/np.sqrt(4*np.pi)*(2*J+1)*np.sqrt(2*k+1)*Wigner_coeff.threej_symbol([[J,k,J],[-M,0,M]])/Wigner_coeff.threej_symbol([[L,k,L],[-L,0,L]])
            pref2 = Wigner_coeff.sixj_symbol([[L,J,S],[J,L,k]])*Wigner_coeff.threej_symbol([[k,3,3],[0,0,0]])
            somma = 0
            for i in range(1,int(conf[1:])+1):
                somma += (-1)**i*Wigner_coeff.threej_symbol([[k,3,3],[0,4-i,i-4]])
            value = pref*pref2*somma*np.sqrt(4*np.pi/(2*k+1))
            coeff.append(value)

    return coeff
