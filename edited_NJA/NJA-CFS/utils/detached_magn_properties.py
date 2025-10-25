from numba import jit, njit
from numba import complex128, boolean, float64, int32, int64
import numba
import numpy as np
import matplotlib.pyplot as plt
from itertools import product, permutations

#======================= DETACHED FUNCTIONS for MAGNETIC PROPERTIES ==============================

@njit(float64(float64))
def fact(number):
    number = int(number)
    if number < 0:
        raise ValueError('Negative number in factorial')
    else:
        factorial = 1.0
        for i in range(1, number + 1, 1):
            factorial *= i
        #print('fact', factorial)
        return factorial

@njit(complex128[:, :](complex128[:, :]))
def from_matrix_to_result_copy(matrix):
    w, v = np.linalg.eig(matrix)
    result = np.zeros((matrix.shape[0] + 1, matrix.shape[0]), dtype=complex128)
    result[0, :] = w
    result[1:, :] = v
    result = result[:, result[0, :].real.argsort()]
    return result

@jit(float64(float64[:,:]))
def sixj_symbol(matrix):
    # computes racah formula for 6j-symbols (p 57 Ch 1 libro Boca)
    # {[a, b, c],[A, B, C]} or {[j1, j2, j3],[j4, j5, j6]}

    # shortcut per calcolare solo quelli che sono != 0
    triads = np.array([[matrix[0,0], matrix[0,1], matrix[0,2]], [matrix[0,0], matrix[1,1], matrix[1,2]],
              [matrix[1,0], matrix[0,1], matrix[1,2]], [matrix[1,0], matrix[1,1], matrix[0,2]]])

    for i in range(len(triads[:,0])):
        if np.sum(matrix[0,:])%1==0:
            pass
        else:
            #print('zero1')
            return 0

        if triads[i,2] >= np.abs(triads[i,0]-triads[i,1]) and triads[i,2] <= triads[i,0]+triads[i,1]:
            pass
        else:
            #print('zero2')
            return 0


    def f(aa, bb, cc):
        # triangular coefficient
        return (fact(aa+bb-cc)*fact(aa-bb+cc)*fact(-aa+bb+cc)/fact(aa+bb+cc+1))**(0.5)

    a = matrix[0,0]
    b = matrix[0,1]
    c = matrix[0,2]
    A = matrix[1,0]
    B = matrix[1,1]
    C = matrix[1,2]

    n_min_list = np.array([0, a+A-c-C, b+B-c-C])
    n_min = max(n_min_list)
    n_max_list = np.array([a+b+A+B+1, a+b-c, A+B-c, a+B-C, A+b-C])
    n_max = min(n_max_list)

    sect1 = (-1)**(a+b+A+B)*f(a,b,c)*f(A,B,c)*f(A,b,C)*f(a,B,C)

    sect2 = 0
    for n in np.arange(n_min,n_max+1,1):
        sect2 += (-1)**n*fact(a+b+A+B+1-n)/(fact(n)*fact(a+b-c-n)*fact(A+B-c-n)*fact(a+B-C-n)*fact(A+b-C-n)*fact(-a-A+c+C+n)*fact(-b-B+c+C+n))

    result = sect1*sect2

    return result

@jit(float64(float64[:,:]))
def threej_symbol(matrix):
    # computes racah formula for 3j-symbols (p 52 Ch 1 libro Boca)
    # ([a , b, c],[A, B, C]) or ([j1, j2, j3],[m1, m2, m3])

    # shortcut per calcolare solo quelli che sono != 0
    for i in range(3):
        if matrix[1,i]>= -matrix[0,i] and matrix[1,i]<= matrix[0,i]:
            pass
        else:
            return 0

    if np.sum(matrix[1,:])==0:
        pass
    else:
        return 0

    if matrix[0,-1] >= np.abs(matrix[0,0]-matrix[0,1]) and matrix[0,-1] <= matrix[0,0]+matrix[0,1]:
        pass
    else:
        return 0

    # if isinstance(np.sum(matrix[0,:]), 'int64')==True:
    if np.sum(matrix[0,:])%1==0:
        if matrix[1,0] == matrix[1,1] and matrix[1,1] == matrix[1,2]:
            if np.sum(matrix[0,:])%2==0:
                pass
            else:
                return 0
        else:
            pass
    else:
        return 0

    a = matrix[0,0]
    b = matrix[0,1]
    c = matrix[0,2]
    A = matrix[1,0]
    B = matrix[1,1]
    C = matrix[1,2]

    n_min_list = np.array([0, -c+b-A, -c+a+B])
    n_min = max(n_min_list)
    n_max_list = np.array([a+b-c, b+B, a-A])
    n_max = min(n_max_list)
    if a-b-C <0:
        factor = (1/(-1))**np.abs(a-b-C)
    else:
        factor = (-1)**(a-b-C)
    sect1 = factor*(fact(a+b-c)*fact(b+c-a)*fact(c+a-b)/fact(a+b+c+1))**(0.5)
    sect2 = (fact(a+A)*fact(a-A)*fact(b+B)*fact(b-B)*fact(c+C)*fact(c-C))**(0.5)
    sect3 = 0
    for n in np.arange(n_min,n_max+1,1):
        sect3 += (-1)**n*(fact(n)*fact(c-b+n+A)*fact(c-a+n-B)*fact(a+b-c-n)*fact(a-n-A)*fact(b-n+B))**(-1)

    result = sect1*sect2*sect3

    # print('dentro',sect1,sect2,sect3,sect1*sect2*sect3)

    return result

@jit(numba.types.UniTuple(complex128[:], 3)(float64,float64,float64,float64,float64,float64,float64,float64,float64[:]))   #capire bene come restituire le tuple
def Zeeman(L,S,J,M,L1,S1,J1,M1,field=np.array([0.,0.,0.])):
    # eq from Boca 2012 p 588

    k=1
    Bx = field[0]
    By = field[1]
    Bz = field[2]

    Bohr = 0.4668604  #conv from BM to cm-1
    ge = 2.0023

    Bq = np.array([-np.sqrt(0.5)*(Bx+1j*By), Bz, np.sqrt(0.5)*(Bx-1j*By)])  # +1, 0, -1

    pre = np.sqrt((2*J+1)*(2*J1+1))

    L1q = k*np.sqrt(L*(L+1)*(2*L+1))*(-1)**(L+S+J+1)*sixj_symbol(np.array([[L, J, S],[J1, L, 1]]))#sixj_symbol(np.array([[J, 1, J1],[L, S, L]]))#sixj_symbol(np.array([[J1, J, 1],[L, L, S]]))

    S1q = ge*np.sqrt(S*(S+1)*(2*S+1))*(-1)**(L+L1+S1*2+J+J1)*(-1)**(L+S+J+1)*sixj_symbol(np.array([[S, J, L],[J1, S, 1]]))#sixj_symbol(np.array([[J, 1, J1],[S, L, S]]))#sixj_symbol(np.array([[J1, J, 1],[S, S, L]]))

    rme = pre*(L1q + S1q)

    int_L1 = np.zeros(3, dtype=complex128)
    int_S1 = np.zeros(3, dtype=complex128)

    integral = 0 + 0 * 1j
    for i, q in enumerate(range(-1, 2, 1)):
        preq = (-1) ** (J - M) * threej_symbol(np.array([[J, 1, J1], [-M, q, M1]]))
        int_L1[i] = preq * pre * L1q
        int_S1[i] = preq * pre * S1q

        integral_Re = (-1) ** q * preq * rme * Bohr * Bq[i].real
        integral_Im = (-1) ** q * preq * rme * Bohr * Bq[i].imag
        integral += integral_Re +1j*integral_Im

    fake_array = np.zeros(3, dtype=complex128)  #this is just because I need to return things of the same type
    fake_array[0] = integral

    return (fake_array, int_L1, int_S1)

@jit(complex128[:,:,:](float64[:,:]))
def mag_moment(basis):
    #costruction of magnetic moment matrix as -kL-geS
    #y component is divided by i (imaginary unit)

    matrix = np.zeros((3, basis.shape[0],basis.shape[0]),dtype=complex128)
    # L_matrix = np.zeros_like(matrix)
    for i in range(basis.shape[0]):
        statei = basis[i]
        Si = statei[0]/2.
        Li = statei[1]
        Ji = statei[2]/2.
        Mi = statei[3]/2.
        seni = statei[4]

        for j in range(0,i+1):
            statej = basis[j]
            Sj = statej[0]/2.
            Lj = statej[1]
            Jj = statej[2]/2.
            Mj = statej[3]/2.
            senj = statej[4]

            if Li==Lj and Si==Sj and seni==senj:
                integral, kL, gS = Zeeman(Li,Si,Ji,Mi,Lj,Sj,Jj,Mj,np.array([0.,0.,0.]))

                # x,y,z  -->  -1,0,+1
                matrix[0,i,j] += (kL[0]+gS[0] - (kL[2]+gS[2]))*1/(np.sqrt(2))
                matrix[1,i,j] += (kL[0]+gS[0] + kL[2]+gS[2])*1j/(np.sqrt(2))
                matrix[2,i,j] += kL[1]+gS[1]

                # L_matrix[0,i,j] += (kL[0] - kL[2])*1/(np.sqrt(2))
                # L_matrix[1,i,j] += (kL[0] + kL[2])*1j/(np.sqrt(2))
                # L_matrix[2,i,j] += kL[1]

                if i!=j:
                    for kk in range(3):
                        matrix[kk,j,i] += np.conj(matrix[kk,i,j])

    matrix = -matrix   #the minus sign is because mu = - kL - 2S
    # print('matrix_mu',matrix[0,...])
    
    return matrix

@jit
def norm(tensor):
    return np.sqrt(np.sum(np.abs(tensor)**2))

@jit
def dfridr(func, x, h, idxi, shape, fargs):

    # print(idxi)

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
            # print('safe exit', a[i, i], a[i - 1, i - 1])
            return result, err

    return result, err

@jit(complex128[:,:](float64[:],float64[:,:],complex128[:,:]))
def add_Zeeman(field_vec, basis, LF_matrix):

    matrix = np.zeros((basis.shape[0],basis.shape[0]),dtype=complex128)
    for i in range(basis.shape[0]):
        statei = basis[i]
        Si = statei[0]/2.
        Li = statei[1]
        Ji = statei[2]/2.
        Mi = statei[3]/2.
        seni = statei[4]

        for j in range(0,i+1):
            statej = basis[j]
            Sj = statej[0]/2.
            Lj = statej[1]
            Jj = statej[2]/2.
            Mj = statej[3]/2.
            senj = statej[4]

            matrix[i,j] += LF_matrix[i,j]
            if i!=j:
                matrix[j,i] += LF_matrix[j,i]

            if Li==Lj and Si==Sj and seni==senj:
                integral, kL, gS = Zeeman(Li,Si,Ji,Mi,Lj,Sj,Jj,Mj,field_vec)
                matrix[i,j] += integral[0]
                if i!=j:
                    matrix[j,i] += np.conj(integral[0])

    return matrix

@jit(float64[:](float64[:],complex128[:,:,:],complex128[:,:],float64[:,:],float64))
def M_vector(field_vec, mu_matrix, LF_matrix, basis, temp):

    kB = 1.380649e-23

    mu = np.zeros((basis.shape[0], 3), dtype=complex128)
    matrix = add_Zeeman(field_vec, basis, LF_matrix)
    result = from_matrix_to_result_copy(matrix)
    E = (result[0,:].real-min(result[0,:].real)) #* 1.9865e-23
    E -= min(E)

    M = np.zeros(3, dtype=float64)

    for kk in range(3):

        den = 0
        num = 0
        for ii in range(len(E)):

            mu_single = np.dot(np.conj(np.ascontiguousarray(result[1:, ii]).T), np.dot(mu_matrix[kk, ...], np.ascontiguousarray(result[1:, ii])))

            if np.abs(mu[ii,kk].imag)<1e-9:
                mu[ii,kk] += mu_single.real
            else:
                print('complex',mu_single)

            num += np.real(mu[ii,kk]*np.exp(-E[ii]/(kB/1.9865e-23*temp)))
            den += np.real(np.exp(-E[ii]/(kB/1.9865e-23*temp)))

        M[kk] = num/den

    return M

@jit('Tuple((float64[:,:],float64[:,:]))(float64[:,:],float64,float64[:,:],complex128[:,:],float64)')
def susceptibility_B_ord1(fields, temp, basis, LF_matrix, delta=0.001):
    # returns the derivative of the function at a point x by Ridders' method of polynomial extrapolation. The value h is input as an estimated initial stepsize.
    # it need not to be small, but rather should be an increment in x over which the function changes substantially. An estimate of the error is also computed.
    # the stepsize is decreased by CON at each iteeration. Max size of tableau is set by NTAB.

    mu0 = 1.25663706212e-06
    muB = 0.4668517532494337

    #print('ord1')
    mu_matrix = mag_moment(basis)  #complex128[:,:,:]
    # print('from ord1: ', mu_matrix)
    chi = np.zeros((fields.shape[0], 3, 3), dtype=float64)
    err = np.zeros_like(chi)
    for i in range(fields.shape[0]):
        for idx in range(3):
            #print('idx',idx)
            chi_comp, err_comp = dfridr(M_vector, fields[i], delta, idx, chi.shape[1:], fargs=(mu_matrix, LF_matrix, basis, temp))
            chi[i,idx] = chi_comp * mu0*muB*1.9865e-23
            err[i,idx] = np.ones(chi_comp.shape)*err_comp * mu0*muB*1.9865e-23

    chi_tensor = np.sum(chi, axis=0)/fields.shape[0]
    err_tensor = np.sum(err, axis=0)/fields.shape[0]

    return (chi_tensor, err_tensor)

def effGval(levels, mu_matrix, v_matrix): 
    #v_matrix = result[1:,:]

    levels = np.array(levels)  
                                                            
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

def calc_torque(B0, T, LF_matrix, basis, plane='xz', step=3.1415/40, figure=False, show_fig=False):

    NA = 6.022e23
    muB = 0.4668517532494337
    Jconv = 1.9865e-23

    def missing_axis(plane_label):

        axes = ['x', 'y', 'z']
        
        missing_axis_label = next(axis for axis in axes if axis not in plane_label)
        missing_axis_index = axes.index(missing_axis_label)
        B_direction = axes.index(plane_label[0])
        
        return missing_axis_label, missing_axis_index, B_direction

    angles = np.arange(0, 2 * np.pi + np.pi/180, step)
    tau_vecs = []
    mu_matrix = mag_moment(basis)

    str_ax, idx_ax, idx_B = missing_axis(plane)
    vec_field = np.array([0.0, 0.0, 0.0])

    # Calculate tau_vec for each B0
    for idx in range(len(B0)):
        vec_field[idx_B] = B0[idx]
        tau_vec = np.zeros((len(angles), 3))
        for i in range(len(angles)):
            rot_vec = rotate_vectors(vec_field, angles[i], str_ax)
            M = M_vector(rot_vec, mu_matrix, LF_matrix, basis, T)
            tau_vec[i, :] = np.cross(M, rot_vec)*NA*muB*Jconv
        tau_vecs.append(tau_vec[:, idx_ax])
    
    if figure:
        fig, ax1 = plt.subplots()
        fig.set_size_inches(7, 4)
        
        line1, = ax1.plot(angles * 180 / np.pi, tau_vecs[0], 'b-', label=f'{B0[0]} T')
        ax1.set_xlabel('Angle / degrees')
        ax1.set_ylabel(r'$\tau$ / Nm/mol', color='b')
        ax1.tick_params(axis='y', labelcolor='b')
        
        axes = [ax1]
        lines = [line1]
        colors = ['r', 'g', 'c', 'm', 'y', 'k']
        for idx in range(1, len(B0)):
            ax = ax1.twinx()
            ax.spines['right'].set_position(('outward', 60 * (idx - 1)))  # Adjust the position
            line, = ax.plot(angles * 180 / np.pi, tau_vecs[idx], color=colors[idx % len(colors)], label=f'{B0[idx]} T')
            ax.set_ylabel(r'$\tau$ / Nm/mol', color=colors[idx % len(colors)])
            ax.tick_params(axis='y', labelcolor=colors[idx % len(colors)])
            axes.append(ax)
            lines.append(line)
        
        labs = [l.get_label() for l in lines]
        ax1.legend(lines, labs, loc='upper right')
        
        ax1.axhline(0, color='black', lw=0.5)
        ax1.set_xlim(0, 360)
        ax1.set_xticks(np.arange(0, 361, 45))
        for i in np.arange(0, 361, 45):
            ax1.axvline(i, color='gray', lw=0.5, alpha=0.5)
        
        ax1.set_title('B direction: ' + plane[0] + ', torque component: ' + str_ax)
        plt.tight_layout()
        if isinstance(figure, str):
            plt.savefig(figure, dpi=600)
        if show_fig:
            plt.show()

    return angles, np.array(tau_vecs)

def susceptibility_field(fields, temp, basis, LF_matrix, delta=0.001):

    kB = 1.380649e-23
    mu0 = 1.25663706212e-06
    muB = 0.4668517532494337

    M_list = np.zeros(fields.shape[0])
    susc_list = np.zeros(fields.shape[0])
    mu = np.zeros((fields.shape[0], basis.shape[0], 3), dtype='complex128')

    mu_matrix = mag_moment(basis)

    for i in range(fields.shape[0]): 

        matrix = add_Zeeman(fields[i], basis, LF_matrix)
        result = from_matrix_to_result_copy(matrix)
        E = (result[0,:].real-min(result[0,:].real)) #* 1.9865e-23   

        den = 0
        num = 0
        for ii in range(len(E)):
            for kk in range(3):
                mu_single = np.dot(np.conj(result[1:,ii]).T, np.dot(mu_matrix[kk,...],result[1:,ii]))
                mu[i,ii,0] += np.copy(mu_single.real)

            num += np.real(mu[i,ii,0]*np.exp(-E[ii]/(kB/1.9865e-23*temp)))
            den += np.real(np.exp(-E[ii]/(kB/1.9865e-23*temp)))
        E_av = num/den/3

        B_inc = fields[i]/np.linalg.norm(fields[i])*delta
        matrix = add_Zeeman(fields[i]+B_inc, basis, LF_matrix)
        result = from_matrix_to_result_copy(matrix)
        E = (result[0,:].real-min(result[0,:].real)) 
        den = 0
        num = 0
        for ii in range(len(E)):
            for kk in range(3):
                mu_single = np.dot(np.conj(result[1:,ii]).T, np.dot(mu_matrix[kk,...],result[1:,ii]))
                mu[i,ii,1] += np.copy(mu_single.real)

            num += np.real(mu[i,ii,1]*np.exp(-E[ii]/(kB/1.9865e-23*temp)))
            den += np.real(np.exp(-E[ii]/(kB/1.9865e-23*temp)))
        E_av_inc = num/den/3

        B_inc = fields[i]/np.linalg.norm(fields[i])*delta
        matrix = add_Zeeman(fields[i]-B_inc, basis, LF_matrix)
        result = from_matrix_to_result_copy(matrix)
        result = from_matrix_to_result(matrix)
        E = (result[0,:].real-min(result[0,:].real)) 
        den = 0
        num = 0
        for ii in range(len(E)):
            for kk in range(3):
                mu_single = np.dot(np.conj(result[1:,ii]).T, np.dot(mu_matrix[kk,...],result[1:,ii]))
                mu[i,ii,2] += np.copy(mu_single.real)

            num += np.real(mu[i,ii,2]*np.exp(-E[ii]/(kB/1.9865e-23*temp)))
            den += np.real(np.exp(-E[ii]/(kB/1.9865e-23*temp)))
        E_av_incm = num/den/3

        M_list[i] = E_av*muB*1.9865e-23
        susc_list[i] = (E_av_inc - E_av_incm)/(2*delta)*mu0*muB*1.9865e-23

    return M_list, susc_list
