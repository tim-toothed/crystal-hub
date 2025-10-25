import numpy as np

#=========================================================================================
#                                    FROM ON_DEV                                         #
#=========================================================================================

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    
def read_AILFT_orca6(filename, conf, method='CASSCF', return_V=False, rotangle_V=False, return_orcamatrix=False):
    # rotate_V is a list of angles to rotate the V matrix
    # rotangle_V = [alpha, beta, gamma] or [q0, q1, q2, q3]
    # euler angles are in scipy convention 'ZYZ'

    from_au = 27.2113834*8065.54477

    file = open(filename).readlines()

    if conf[0] == 'f':
        collected = []
        check_name = []
        passed = False
        F2 = 0
        F4 = 0
        F6 = 0
        zeta = 0
        for i, line in enumerate(file):
            if "AILFT MATRIX ELEMENTS ("+method+")" in line:
                for j in range(4, 11):
                    collected.append(file[i+j][10:])
                    check_name.append(file[i+j][:10])
                    passed = True
            if passed and not F2 and 'F2' in line:
                splitline = line.split(' ')
                numbers = [float(num) for num in splitline if is_float(num)]
                F2 = numbers[-1]
            if passed and not F4 and 'F4' in line:
                splitline = line.split(' ')
                numbers = [float(num) for num in splitline if is_float(num)]
                F4 = numbers[-1]
            if passed and not F6 and 'F6' in line:
                splitline = line.split(' ')
                numbers = [float(num) for num in splitline if is_float(num)]
                F6 = numbers[-1]
            if passed and not zeta and 'SOC constant zeta' in line:
                splitline = line.split(' ')
                numbers = [float(num) for num in splitline if is_float(num)]
                zeta = numbers[-1]
                passed = False  #since in orca output zeta is usually after F2

        names = ['f0', 'f+1', 'f-1', 'f+2', 'f-2', 'f+3', 'f-3']
        for i,item in enumerate(check_name):
            if names[i] not in item:
                raise ValueError(f"Error in reading the file: {names[i]} not found in {item} while reading ORCA file {filename}")
        matrix = []
        for i in range(len(collected)):
            splitline = collected[i].split(' ')
            row = []
            for j in range(len(splitline)):
                if splitline[j] != '' and is_float(splitline[j]):
                    row.append(float(splitline[j]))
            matrix.append(row)
        matrix = np.array(matrix) #order 0 1 -1 2 -2 3 -3 

        if rotangle_V:
            if len(rotangle_V) > 3:
                matrix = rotate_V_quat(matrix, 3, rotangle_V)
            else:
                matrix = rotate_V(matrix, 3, *rotangle_V)

        matrix = adjust_phase_matrix(matrix)

        if return_orcamatrix:
            return matrix

        ### change to order 0 -1 1 -2 2 -3 3 
        mask = [[(0,0), (2,0), (1,0), (4,0), (3,0), (6,0), (5,0)],
                [(2,0), (2,2), (1,2), (4,2), (3,2), (6,2), (5,2)],
                [(1,0), (1,2), (1,1), (4,1), (3,1), (6,1), (5,1)],
                [(4,0), (4,2), (4,1), (4,4), (3,4), (6,4), (5,4)],
                [(3,0), (3,2), (3,1), (3,4), (3,3), (6,3), (5,3)],
                [(6,0), (6,2), (6,1), (6,4), (6,3), (6,6), (5,6)],
                [(5,0), (5,2), (5,1), (5,4), (5,3), (5,6), (5,5)]]
        mask = np.array(mask)
        row_indices = mask[:,:,0]
        col_indices = mask[:,:,1]
        matrix = matrix[row_indices, col_indices]

        dic = {}

        dic_V = {
            '11':matrix[0,0]*from_au,
            '21':matrix[1,0]*from_au, '22':matrix[1,1]*from_au,
            '31':matrix[2,0]*from_au,'32':matrix[2,1]*from_au,'33':matrix[2,2]*from_au,
            '41':matrix[3,0]*from_au,'42':matrix[3,1]*from_au,'43':matrix[3,2]*from_au,'44':matrix[3,3]*from_au,
            '51':matrix[4,0]*from_au,'52':matrix[4,1]*from_au,'53':matrix[4,2]*from_au,'54':matrix[4,3]*from_au,'55':matrix[4,4]*from_au,
            '61':matrix[5,0]*from_au,'62':matrix[5,1]*from_au,'63':matrix[5,2]*from_au,'64':matrix[5,3]*from_au,'65':matrix[5,4]*from_au,'66':matrix[5,5]*from_au,
            '71':matrix[6,0]*from_au,'72':matrix[6,1]*from_au,'73':matrix[6,2]*from_au,'74':matrix[6,3]*from_au,'75':matrix[6,4]*from_au,'76':matrix[6,5]*from_au,'77':matrix[6,6]*from_au
        }
        
        dic_Bkq = from_Vint_to_Bkq(dic_V, conf)
        dic['dic_bkq'] = dic_Bkq
        dic['F2'] = F2
        dic['F4'] = F4
        dic['F6'] = F6
        dic['zeta'] = zeta


    else:

        collected = []
        check_name = []
        passed = False
        F2 = 0
        F4 = 0
        zeta = 0
        for i, line in enumerate(file):
            if "AILFT MATRIX ELEMENTS ("+method+")" in line:
                for j in range(5, 10):
                    collected.append(file[i+j][10:])
                    check_name.append(file[i+j][:10])
                    passed = True
            if passed and not F2 and 'F2' in line:
                splitline = line.split(' ')
                numbers = [float(num) for num in splitline if is_float(num)]
                F2 = numbers[-1]
            if passed and not F4 and 'F4' in line:
                splitline = line.split(' ')
                numbers = [float(num) for num in splitline if is_float(num)]
                F4 = numbers[-1]
            if passed and not zeta and 'SOC constant zeta' in line:
                splitline = line.split(' ')
                numbers = [float(num) for num in splitline if is_float(num)]
                zeta = numbers[-1]
                passed = False  #since in orca output zeta is usually after F2

        names = ['dz2', 'dxz', 'dyz', 'dx2-y2', 'dxy']
                  
        for i,item in enumerate(check_name):
            if names[i] not in item:
                raise ValueError(f"Error in reading the file: {names[i]} not found in {item} while reading ORCA file {filename}")
        matrix = []
        for i in range(len(collected)):
            splitline = collected[i].split(' ')
            row = []
            for j in range(len(splitline)):
                if splitline[j] != '' and is_float(splitline[j]):
                    row.append(float(splitline[j]))
            matrix.append(row)
        matrix = np.array(matrix)

        if rotangle_V:
            if len(rotangle_V) > 3:
                matrix = rotate_V_quat(matrix, 2, rotangle_V)
            else:
                matrix = rotate_V(matrix, 2, *rotangle_V)

        if return_orcamatrix:
            return matrix

        ### change to order  z2, yz, xz, xy, x2-y2 
        mask = [[(0,0), (2,0), (1,0), (4,0), (3,0)],
                [(2,0), (2,2), (1,2), (4,2), (3,2)],
                [(1,0), (1,2), (1,1), (4,1), (3,1)],
                [(4,0), (4,2), (4,1), (4,4), (3,4)],
                [(3,0), (3,2), (3,1), (3,4), (3,3)]]
        mask = np.array(mask)
        row_indices = mask[:,:,0]
        col_indices = mask[:,:,1]
        matrix = matrix[row_indices, col_indices]

        dic = {}

        dic_V = {
            '11':matrix[0,0]*from_au,
            '21':matrix[1,0]*from_au, '22':matrix[1,1]*from_au,
            '31':matrix[2,0]*from_au,'32':matrix[2,1]*from_au,'33':matrix[2,2]*from_au,
            '41':matrix[3,0]*from_au,'42':matrix[3,1]*from_au,'43':matrix[3,2]*from_au,'44':matrix[3,3]*from_au,
            '51':matrix[4,0]*from_au,'52':matrix[4,1]*from_au,'53':matrix[4,2]*from_au,'54':matrix[4,3]*from_au,'55':matrix[4,4]*from_au,
        }
        
        dic_Bkq = from_Vint_to_Bkq(dic_V, conf)
        dic['dic_bkq'] = dic_Bkq
        dic['F2'] = F2
        dic['F4'] = F4
        dic['zeta'] = zeta

    if return_V:
        return dic_V
    else:
        return dic

def adjust_phase_matrix(matrix):
    #adjust the phase of the LF matrix from ORCA
    for i in range(7):  # Python uses 0-based indexing
        phasei = 1.0
        if i == 5 or i == 6:  # Adjusted for 0-based indexing
            phasei = -1.0
        for j in range(7):
            phasej = 1.0
            if j == 5 or j == 6:  # Adjusted for 0-based indexing
                phasej = -1.0
            matrix[i, j] *= phasei * phasej
    return matrix

def rotate_V(V, l, alpha=0.0, beta=0.0, gamma=0.0, real=True):
    #in rad
    #V has to be in the order 0, 1, -1, 2, -2, 3, -3 (like directly from orca)
    V_rot = np.zeros_like(V, dtype='complex128')
    M_num = [0, 1, -1, 2, -2, 3, -3]
    M1_num = [0, 1, -1, 2, -2, 3, -3]
    U = U_complex2real(l)
    D = np.zeros_like(V, dtype='complex128')
    for im1, m1 in enumerate(np.arange(-l, l+1, 1)):
        for jm, m in enumerate(np.arange(-l, l+1, 1)):
            D[im1, jm] = Wigner_coeff.Wigner_Dmatrix(l, m1, m, alpha, beta, gamma)
    if real:
        R = np.conj(U).T@D@U   
    else:
        R = D.copy()
    mask = []
    perm = [3,4,2,5,1,6,0]
    for i,ii in enumerate(perm):
        row = []
        for j,jj in enumerate(perm):
            row.append((ii,jj))
        mask.append(row)
    mask = np.array(mask)
    row_indices = mask[:,:,0]
    col_indices = mask[:,:,1]
    R = R[row_indices, col_indices]
    for iM,M in enumerate(M_num):
        for jM,M1 in enumerate(M1_num):
            for im1,m1 in enumerate(np.arange(-l, l+1, 1)):
                for jm2,m2 in enumerate(np.arange(-l, l+1, 1)):
                    V_rot[iM,jM] += R[iM,im1]*V[im1,jm2]*R[jM,jm2]
    V_rot = V_rot.real
    return V_rot

def rotate_V_quat(V, l, q=[1.0,0.0,0.0,0.0], real=True):
    #V has to be in the order 0, 1, -1, 2, -2, 3, -3 (like directly from orca)
    V_rot = np.zeros_like(V, dtype='complex128')
    M_num = [0, 1, -1, 2, -2, 3, -3]
    M1_num = [0, 1, -1, 2, -2, 3, -3]
    dict, coeff = read_DWigner_quat()
    U = U_complex2real(l)
    D = Wigner_coeff.Wigner_Dmatrix_quat_complete(l, q, dict = dict, coeff = coeff)   #order: -3,-2,-1,0,1,2,3
    if real:
        R = np.conj(U).T@D@U
    else:
        R = D.copy()
    mask = []
    perm = [3,4,2,5,1,6,0]
    for i,ii in enumerate(perm):
        row = []
        for j,jj in enumerate(perm):
            row.append((ii,jj))
        mask.append(row)
    mask = np.array(mask)
    row_indices = mask[:,:,0]
    col_indices = mask[:,:,1]
    R = R[row_indices, col_indices]
    for iM,M in enumerate(M_num):
        for jM,M1 in enumerate(M1_num):
            for im1,m1 in enumerate(np.arange(-l, l+1, 1)):
                for jm2,m2 in enumerate(np.arange(-l, l+1, 1)):
                    V_rot[iM,jM] += R[iM,im1]*V[im1,jm2]*np.conj(R[jM,jm2])
    V_rot = V_rot.real
    return V_rot

def rotate_dicV(dic, l, rotangle_V=False, real=True, return_orcaV = False):
    #ruota V a partire dal dizionario nell'ordine 0, -1, 1, -2, 2, -3, 3
    #quindi prima lo converte in matrice 0, 1, -1, 2, -2, 3, -3
    #poi lo ruota
    #poi lo converte in dizionario 0, -1, 1, -2, 2, -3, 3

    def permutation(V):
        perm = [0,2,1,4,3,6,5]
        mask = []
        for i,ii in enumerate(perm):
            row = []
            for j,jj in enumerate(perm):
                row.append((ii,jj))
            mask.append(row)
        mask = np.array(mask)
        row_indices = mask[:,:,0]
        col_indices = mask[:,:,1]
        V = V[row_indices, col_indices]
        return V
    
    #-----------------------------------------------------------------

    V = np.zeros((2*l+1, 2*l+1))
    for i in range(2*l+1):
        for j in range(0,i+1):
            V[i,j] = dic[str(i+1)+str(j+1)]
            V[j,i] = dic[str(i+1)+str(j+1)]
    
    #da [0, -1, 1, -2, 2, -3, 3] a [0, 1, -1, 2, -2, 3, -3]
    V = permutation(V)
    if return_orcaV:
        return V/27.2113834/8065.54477

    #ruoto
    if rotangle_V is not False:
        if len(rotangle_V) > 3:
            V_rot = rotate_V_quat(V, l, rotangle_V, real)
        else:
            V_rot = rotate_V(V, l, *rotangle_V, real)

    #da [0, 1, -1, 2, -2, 3, -3] a [0, -1, 1, -2, 2, -3, 3]
    V_rot = permutation(V_rot)

    dic_rot = {}
    for i in range(2*l+1):
        for j in range(0,i+1):
            dic_rot[str(i+1)+str(j+1)] = V_rot[i,j]

    return dic_rot
    
def U_complex2real(l):
    # order: ... -3, -2, -1, 0, 1, 2, 3 ...
    dim = 2*l + 1
    U = 1j * np.zeros((dim, dim))
    mvalues = range(l, -l - 1, -1)
    for i, m in enumerate(mvalues):
        if m < 0:
            U[i, i] = (-1)**m / np.sqrt(2)
            U[dim - i - 1, i] = 1 / np.sqrt(2)
        elif m == 0:
            U[i, i] = 1.0
        else:
            U[i, i] = 1j / np.sqrt(2)
            U[dim - i - 1, i] = -(-1)**m * 1j / np.sqrt(2)
    return U

def read_xyz(filename):
    #the first two lines must be the number of atoms and a comment line
    #in the subsequent lines the first column must be the atom name

    filecoord = []
    for i,line in enumerate(open(filename).readlines()):
        if i>1:
            splitline = line.split()
            num = [float(num) for num in splitline if is_float(num)]
            filecoord.append([splitline[0]]+num)

    filecoord = np.array(filecoord)
    return filecoord

def Rzyz_active(alpha=0,beta=0,gamma=0):
    #proper-euler angles
    #active rotation

    ca = np.cos(alpha)
    cb = np.cos(beta)
    cc = np.cos(gamma)
    sa = np.sin(alpha)
    sb = np.sin(beta)
    sc = np.sin(gamma)
    
    R1 = np.array([ca*cb*cc-(sa*sc), -ca*cb*sc-(sa*cc), ca*sb])
    R2 = np.array([ca*sc+(sa*cb*cc), -sa*cb*sc+(ca*cc), sa*sb])
    R3 = np.array([-sb*cc, sb*sc, cb])
    R = np.stack((R1, R2, R3), axis=0)
    
    return R 
