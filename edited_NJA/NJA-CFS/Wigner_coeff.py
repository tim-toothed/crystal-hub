import numpy as np

class Wigner_coeff():
    """
    This class is a wrapper around the functions to compute angular momentum coupling coefficients
    and Wigner D matrices for spherical harmonics rotations.

    References:
    1. Boca, R., "Theoretical Foundations of Molecular Magnetism" (1999)
    """

    def fact(number):
        """
        Calculate the factorial of a given number.

        This function takes an integer as input and calculates its factorial. 
        The factorial of a number is the product of all positive integers less than or equal to the number.
        If a negative number is provided, the function prints an error message and exits.

        Parameters:
        number (int): The number to calculate the factorial of.

        Returns:
        factorial (int): The factorial of the input number.
        """
        number = int(number)
        if number < 0:
            print('negative number in factorial')
            exit()
        else:
            factorial = 1
            for i in range(1,number+1,1):
                factorial *= i
            return int(factorial)

    def threej_symbol(matrix):
        # computes racah formula for 3j-symbols (p 52 Ch 1 from Boca)
        # ([a , b, c],[A, B, C]) or ([j1, j2, j3],[m1, m2, m3])
        """
        Compute the 3j-symbol using the Racah formula.

        This function takes a 2D array as input, representing a 3j-symbol in the form ([j1, j2, j3],[m1, m2, m3]).
        It then calculates the 3j-symbol using the Racah formula, which is used in quantum mechanics to calculate the coefficients in the transformation of a couple of angular momenta.

        Parameters:
        matrix (numpy.ndarray): A 2D array representing a 3j-symbol.

        Returns:
        result (float): The calculated 3j-symbol.
        """

        matrix = np.array(matrix)

        # print(matrix)

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

        n_min_list = [0, -c+b-A, -c+a+B]
        n_min = max(n_min_list)
        n_max_list = [a+b-c, b+B, a-A]
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

        return result

    def sixj_symbol(matrix):
        # computes racah formula for 6j-symbols (p 57 Ch 1 libro Boca)
        # {[a, b, c],[A, B, C]} or {[j1, j2, j3],[j4, j5, j6]}
        """
        Compute the 6j-symbol using the Racah formula.

        This function takes a 2D array as input, representing a 6j-symbol in the form ([j1, j2, j3],[j4, j5, j6]).
        It then calculates the 6j-symbol using the Racah formula, which is used in quantum mechanics to calculate the coefficients in the transformation of a couple of angular momenta.

        Parameters:
        matrix (numpy.ndarray): A 2D array representing a 6j-symbol.

        Returns:
        result (float): The calculated 6j-symbol.
        """

        matrix = np.array(matrix)
        #print(matrix)

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
            res = (fact(aa+bb-cc)*fact(aa-bb+cc)*fact(-aa+bb+cc)/fact(aa+bb+cc+1))**(0.5)
            if res.imag!=0:
                return 0
            else:
                return res

        a = matrix[0,0]
        b = matrix[0,1]
        c = matrix[0,2]
        A = matrix[1,0]
        B = matrix[1,1]
        C = matrix[1,2]

        n_min_list = [0, a+A-c-C, b+B-c-C]
        n_min = max(n_min_list)
        n_max_list = [a+b+A+B+1, a+b-c, A+B-c, a+B-C, A+b-C]
        n_max = min(n_max_list)

        sect1 = (-1)**(a+b+A+B)*f(a,b,c)*f(A,B,c)*f(A,b,C)*f(a,B,C)

        #print(sect1)

        sect2 = 0
        for n in np.arange(n_min,n_max+1,1):
            sect2 += (-1)**n*fact(a+b+A+B+1-n)/(fact(n)*fact(a+b-c-n)*fact(A+B-c-n)*fact(a+B-C-n)*fact(A+b-C-n)*fact(-a-A+c+C+n)*fact(-b-B+c+C+n))
            #print(sect2)

        result = sect1*sect2

        return result

    def Wigner_Dmatrix(l, m1, m, alpha, beta, gamma):
        #calcola un elemento della matrice di wigner D^l_m1m(alpha,beta,gamma)
        #gli angoli devono essere dati in radianti
        """
        Compute an element of the Wigner D-matrix (D_{m1,m}^{l}).

        This function calculates an element of the Wigner D-matrix, which is used in quantum mechanics to describe rotations of quantum states.
        The Wigner D-matrix is defined in terms of Euler angles.

        Parameters:
        l (int): The one-electron orbital quantum number.
        m1 (int): The initial magnetic quantum number.
        m (int): The final magnetic quantum number.
        alpha (float): The first Euler angle, in radians.
        beta (float): The second Euler angle, in radians.
        gamma (float): The third Euler angle, in radians.

        Returns:
        result (complex): The calculated element of the Wigner D-matrix.
        """

        d = np.sqrt(fact(l+m1)*fact(l-m1)*fact(l+m)*fact(l-m))
        somma = 0
        smin = max(0, m-m1)
        smax = min(l+m,l-m1)
        for s in range(smin,smax+1):
            somma += (-1)**(m1-m+s)*np.cos(beta/2)**(2*l+m-m1-2*s)*np.sin(beta/2)**(m1-m+2*s)/(fact(l+m-s)*fact(s)*fact(m1-m+s)*fact(l-m1-s))
        d *= somma
        return np.exp(-1j*m1*alpha)*d*np.exp(-1j*m*gamma)

    def Wigner_Dmatrix_quat_complete(l, R, bin=1e-8, dict = None, coeff=None):
        #R = R1 1 + Rx x + Ry y + Rz z
        #equations from reorientational correlation functions, quaternions and wigner rotation matrices
        # print(R)
        """
        Compute the Wigner D-matrix using quaternions: R = R1 1 + Rx x + Ry y + Rz z.

        This function calculates the Wigner D-matrix using quaternions, which are a more efficient way to represent rotations than Euler angles.
        The Wigner D-matrix is used in quantum mechanics to describe rotations of quantum states.
        The quaternion in input does not need to be normalized since it is normalized before the calculation.

        Parameters:
        l (int): The one-electron orbital quantum number.
        R (numpy.ndarray): The quaternion representing the rotation.
        bin (float, optional): The tolerance size for the calculation. Default is 1e-8.
        dict (dict, optional): dictionary from tables. Default is None.
        coeff (numpy.ndarray, optional): The coefficients for the calculation from table. Default is None.

        Returns:
        D (numpy.ndarray): The calculated Wigner D-matrix as a complex numpy array.
        """
        R /= np.linalg.norm(R)

        A = R[0] -1j*R[3]
        B = -R[2]-1j*R[1]

        if np.abs(A)<bin:
            A = 0
        if np.abs(B)<bin:
            B = 0

        Z = R[0]**2 - R[1]**2 - R[2]**2 + R[3]**2
        Ac = np.conj(A)
        Bc = np.conj(B)
        D = np.zeros((2*l+1, 2*l+1), dtype='complex128')
        if l==2:
            D[4,4] = A**4                           # D[2,2] =
            D[0,0] = np.conj(D[4,4])              #D[-2,-2]
            D[4,3] = 2*A**3*B                       #D[2,1] =
            D[0,1] = (-1)*np.conj(D[4,3])         #D[-2,-1]
            D[4,2] = np.sqrt(6)*A**2*B**2           #D[2,0] =
            D[0,2] = np.conj(D[4,2])               #D[-2,0] =
            D[3,4] = -2*A**3*Bc                     #D[1,2] =
            D[1,0] = -1*np.conj(D[3,4])           #D[-1,-2]
            D[3,3] = A**2*(2*Z-1)                   #D[1,1] =
            D[1,1] = np.conj(D[3,3])              #D[-1,-1]
            D[3,2] = np.sqrt(6)*A*B*Z               #D[1,0] =
            D[1,2] = -1*np.conj(D[3,2])            #D[-1,0] =
            D[2,4] = np.sqrt(6)*A**2*Bc**2          #D[0,2] =
            D[2,0] = np.conj(D[2,4])               #D[0,-2] =
            D[2,3] = -np.sqrt(6)*A*Bc*Z             #D[0,1] =
            D[2,1] = -1*np.conj(D[2,3])            #D[0,-1] =
            D[2,2] = 1/2 * (3*Z**2-1)               #D[0,0] =
            D[1,4] = -2*A*Bc**3                    #D[-1,2] =
            D[3,0] = -np.conj(D[1,4])             #D[1,-2] =
            D[1,3] = Bc**2*(2*Z+1)                 #D[-1,1] =
            D[3,1] = np.conj(D[1,3])              #D[1,-1] =
            D[0,4] = Bc**4                         #D[-2,2] =
            D[4,0] = np.conj(D[0,4])              #D[2,-2] =
            D[0,3] = -2*Ac*Bc**3                   #D[-2,1] =
            D[4,1] = -np.conj(D[0,3])             #D[2,-1] =
        if l==3:
            D[6,6] = A**6                                       #D[3,3] = A
            D[0,0] = Ac**6                                      #D[-3,-3] =
            D[6,5] = np.sqrt(6)*B*A**5                          #D[3,2] = n
            D[0,1] = -np.sqrt(6)*Bc*Ac**5                         #D[-3,-2] =
            D[6,4] = np.sqrt(15)*B**2*A**4                      #D[3,1] = n
            D[0,2] = np.sqrt(15)*Bc**2*Ac**4                  #D[-3,-1] =
            D[6,3] = 2*np.sqrt(5)*B**3*A**3                     #D[3,0] = 2
            D[0,3] = -2*np.sqrt(5)*Bc**3*Ac**3                #D[-3,0] =
            D[5,6] = -np.sqrt(6)*A**5*Bc                        #D[2,3] = -
            D[1,0] = np.sqrt(6)*Ac**5*B                       #D[-2,-3] =
            D[5,5] = A**4*(3*Z-2)                               #D[2,2] = A
            D[1,1] = Ac**4*(3*Z-2)                            #D[-2,-2] =
            D[5,4] = 1/2 * np.sqrt(10)*B*A**3*(3*Z-1)           #D[2,1] = 1
            D[1,2] = -1/2*np.sqrt(10)*Bc*Ac**3*(3*Z-1)        #D[-2,-1] =
            D[5,3] = np.sqrt(30)*B**2*A**2*Z                    #D[2,0] = n
            D[1,3] = np.sqrt(30)*Bc**2*Ac**2*Z                 #D[-2,0] =
            D[4,6] = np.sqrt(15)*A**4*Bc**2                     #D[1,3] = n
            D[2,0] = np.sqrt(15)*Ac**4*B**2                   #D[-1,-3] =
            D[4,5] = 1/2 * np.sqrt(10)*A**3*Bc*(1-3*Z)          #D[1,2] = 1
            D[2,1] = -1/2*np.sqrt(10)*Ac**3*B*(1-3*Z)         #D[-1,-2] =
            D[4,4] = 1/4*A**2*(15*Z**2-10*Z-1)                  #D[1,1] = 1
            D[2,2] = 1/4*Ac**2*(15*Z**2-10*Z-1)               #D[-1,-1] =
            D[4,3] = 1/2 * np.sqrt(3)*A*B*(5*Z**2-1)            #D[1,0] = 1
            D[2,3] = -1/2 * np.sqrt(3)*Ac*Bc*(5*Z**2-1)        #D[-1,0] =
            D[3,6] = -2*np.sqrt(5)*A**3*Bc**3                   #D[0,3] = -
            D[3,0] = 2*np.sqrt(5)*Ac**3*B**3                   #D[0,-3] =
            D[3,5] = np.sqrt(30)*A**2*Bc**2*Z                   #D[0,2] = n
            D[3,1] = np.sqrt(30)*Ac**2*B**2*Z                  #D[0,-2] =
            D[3,4] = 1/2 * np.sqrt(3)*A*Bc*(1-5*Z**2)          #D[0,1] = -
            D[3,2] = -1/2*np.sqrt(3)*Ac*B*(1-5*Z**2)            #D[0,-1] =
            D[3,3] = 1/2*(5*Z**3-3*Z)                           #D[0,0] = 1
            D[2,6] = np.sqrt(15)*A**2*Bc**4                    #D[-1,3] =
            D[4,0] = np.sqrt(15)*Ac**2*B**4                    #D[1,-3] =
            D[2,5] = -1/2 * np.sqrt(10)*A*Bc**3*(1+3*Z)        #D[-1,2] =
            D[4,1] = 1/2*np.sqrt(10)*Ac*B**3*(1+3*Z)           #D[1,-2] =
            D[2,4] = 1/4 * Bc**2*(15*Z**2+10*Z-1)              #D[-1,1] =
            D[4,2] = 1/4 * B**2*(15*Z**2+10*Z-1)               #D[1,-1] =
            D[1,6] = -np.sqrt(6)*A*Bc**5                       #D[-2,3] =
            D[5,0] = np.sqrt(6)*Ac*B**5                        #D[2,-3] =
            D[1,5] = Bc**4*(3*Z+2)                             #D[-2,2] =
            D[5,1] = B**4*(3*Z+2)                              #D[2,-2] =
            D[1,4] = -1/2 * np.sqrt(10)*Ac*Bc**3*(3*Z+1)       #D[-2,1] =
            D[5,2] = 1/2*np.sqrt(10)*A*B**3*(3*Z+1)            #D[2,-1] =
            D[0,6] = Bc**6                                     #D[-3,3] =
            D[6,0] = B**6                                      #D[3,-3] =
            D[0,5] = -np.sqrt(6)*Ac*Bc**5                      #D[-3,2] =
            D[6,1] = np.sqrt(6)*A*B**5                         #D[3,-2] =
            D[0,4] = np.sqrt(15)*Bc**4*Ac**2                   #D[-3,1] =
            D[6,2] = np.sqrt(15)*B**4*A**2                     #D[3,-1] =
        elif (l==4 or l==6) and dict is not None and coeff is not None:
            part = dict[l]
            for key,value in part.items():
                idx = [int(ii) for ii in key.split(':')]
                i = np.abs(idx[0]+l)
                j = np.abs(idx[1]+l)
                D[i,j] = eval(value)
                D[i,j] *= coeff[l][i,j]

            D_rep = np.zeros_like(D, dtype='complex128')
            for i,ii in enumerate(range(l,-1,-1)):
                for j,jj in enumerate(range(l,-l-1,-1)):
                    D_rep[i,j] = D[i,j] 
                    D_rep[-i-1,-j-1] = (-1)**(np.abs(i-j))*np.conj(D[i,j]) 
            D=D_rep

        return D
