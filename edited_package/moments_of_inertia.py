"""
Geometric symmetry detection and moment of inertia analysis for ligand coordination.

Identifies rotational symmetry axes (n-fold, n=2,3,4,5) in ligand geometries using
Continuous Shape Measures (CSM). Falls back to moment of inertia tensor when no
clear symmetry exists. Used to automatically orient crystal field calculations.

Key algorithms:
- CSM-based symmetry detection (Pinsky & Avnir, Inorg. Chem. 37, 5575, 1998)
- Moment of inertia principal axes
- Rodrigues rotation formula with numba acceleration
"""

import numpy as np
from numba import njit
from scipy.optimize import minimize

# Outer function
def findZaxis(atoms):
    """
    Automatically detect primary symmetry axis in ligand geometry.
    
    Tests for n-fold rotation axes (n=4,3,5,2) using Continuous
    Shape Measures threshold < 1.
    If no symmetry found, uses moment of inertia outlier axis.
    
    Args:
        atoms: Array of ligand position vectors (Nx3)
    
    Returns:
        tuple: (z_axis, y_axis) - orthogonal basis vectors
    
    Note:
        Priority order: 4-fold > 3-fold > 5-fold > 2-fold > inertia
    """

    RotA4, f4 = findZaxis_SOM_rotation(atoms, np.pi/2)  ## four-fold rotation
    # print('############ FOUR FOLD ROTATION CSM =',f4,'\n')
    if f4 < 1:
        yax = np.cross(RotA4, atoms[0])
        yax /= np.linalg.norm(yax)
        print('\tFound a near-4-fold axis...  CSM=', f4)
        return RotA4, yax
    else:
        RotA3, f3 = findZaxis_SOM_rotation(atoms, np.pi/3*2)  ## three-fold rotation
        if f3 < 1:
            yax = np.cross(RotA3, atoms[0])
            yax /= np.linalg.norm(yax)
            print('\tFound a near-3-fold axis...  CSM=', f3)
            return RotA3, yax
        else:
            RotA5, f5 = findZaxis_SOM_rotation(atoms, np.pi/5*2)  ## five-fold rotation
            if f5 < 1:
                yax = np.cross(RotA5, atoms[0])
                yax /= np.linalg.norm(yax)
                print('\tFound a near-5-fold axis...  CSM=', f5)
                return RotA3, yax
            else:
                RotA2, f2 = findZaxis_SOM_rotation(atoms, np.pi)  ## two-fold rotation
                if f2 < 1:
                    yax = np.cross(RotA2, atoms[0])
                    yax /= np.linalg.norm(yax)
                    print('\tFound a near-2-fold axis...  CSM=', f2)
                    return RotA2, yax
                else:
                    print('\tUsing moment of intertia tensor to estimate z axis...')
                    return selectZaxisMI(atoms) ## Select using moment of intertia



# Inner function
def MomIntertia(atoms):
    """
    Calculate moment of inertia tensor.
    
    Computes 3x3 tensor I with elements:
        I_xx = Σ(y² + z²), I_xy = -Σ(xy), etc.
    
    Args:
        atoms: Atomic position vectors (Nx3)
    
    Returns:
        np.ndarray: 3x3 symmetric inertia tensor
    """
    II = np.zeros((3,3))
    for at in atoms:
        x,y,z = at
        II[0,0] += y**2 + z**2
        II[1,1] += x**2 + z**2
        II[2,2] += x**2 + y**2
        II[0,1] += -x*y
        II[0,2] += -x*z
        II[1,2] += -y*z
        II[1,0] += -x*y
        II[2,0] += -x*z
        II[2,1] += -y*z
    return II

# Inner function
def selectZaxisMI(atoms):
    """
    Select z-axis from moment of inertia eigenanalysis.
    
    Chooses principal axis with most distinct eigenvalue as z.
    Y-axis is most similar eigenvalue (disk-like distribution).
    
    Args:
        atoms: Position vectors (Nx3)
    
    Returns:
        tuple: (z_axis, y_axis) eigenvectors
    """
    II = MomIntertia(atoms)
    evals, evecs = np.linalg.eig(II)
    # Find the "outlier" eigenvalue
    evd = [np.abs(evals[0] - evals[1]) + np.abs(evals[0] - evals[2]),
           np.abs(evals[1] - evals[0]) + np.abs(evals[1] - evals[2]),
           np.abs(evals[2] - evals[0]) + np.abs(evals[2] - evals[1])]
    rotAxis = evecs.T[np.argmax(evd)]
    yaxis = evecs.T[np.argmin(evd)]
    ## Ensure that the y axis is orthogonal...
    ## yaxis -= np.dot(yaxis,rotAxis)/np.linalg.norm(rotAxis)
    return rotAxis, yaxis

# Inner function
def ContinuousShapeMeasure(shape1,shape2):
    """
    Calculate CSM between two point sets.
    
    CSM = Σ_i min_j |r_i - s_j|²
    Measures similarity: 0 = identical, larger = more different
    
    Args:
        shape1, shape2: Point sets (Nx3 arrays)
    
    Returns:
        float: Shape measure value
    """
    CSM = 0
    for r1 in np.array(shape1):
        dist = np.sum((np.array(shape2) - r1)**2, axis=1)
        CSM += np.min(dist)
    return CSM

# Inner function
@njit
def anglesToVector(theta,phi):
    """
    Convert spherical angles to Cartesian unit vector.
    
    Args:
        theta: Polar angle (0 to π)
        phi: Azimuthal angle (0 to 2π)
    
    Returns:
        tuple: (x, y, z) unit vector
    """
    return np.sin(theta)*np.cos(phi),  np.sin(theta)*np.sin(phi),  np.cos(theta)

# Inner function
@njit
def rotationMatrix(theta,phi, angle):
    """
    Rodrigues rotation matrix about axis (θ,φ) by angle.
    
    Numba-accelerated for optimization loops.
    
    Args:
        theta, phi: Rotation axis in spherical coords
        angle: Rotation angle (radians)
    
    Returns:
        np.ndarray: 3x3 rotation matrix
    """
    u, v, w = anglesToVector(theta,phi)
    rotmatrix = np.zeros((3,3))
    rotmatrix[0,0] = (u**2 +(v**2 + w**2)*np.cos(angle))
    rotmatrix[0,1] = (u*v*(1- np.cos(angle)) - w*np.sin(angle))
    rotmatrix[0,2] = (u*w*(1- np.cos(angle)) + v*np.sin(angle))
    rotmatrix[1,0] = (u*v*(1- np.cos(angle)) + w*np.sin(angle))
    rotmatrix[1,1] = (v**2 +(u**2 + w**2)*np.cos(angle))
    rotmatrix[1,2] = (v*w*(1- np.cos(angle)) - u*np.sin(angle))
    rotmatrix[2,0] = (u*w*(1- np.cos(angle)) - v*np.sin(angle))
    rotmatrix[2,1] = (v*w*(1- np.cos(angle)) + u*np.sin(angle))
    rotmatrix[2,2] = (w**2 +(v**2 + u**2)*np.cos(angle))
    return rotmatrix

# Inner function
def rotateArbAxis(atoms, theta,phi, angle):
    """
    Apply rotation to point set.
    
    Args:
        atoms: Position vectors (Nx3)
        theta, phi: Axis direction
        angle: Rotation angle
    
    Returns:
        list: Rotated positions
    """
    rotmat = rotationMatrix(theta,phi, angle)
    newat = []
    for at in atoms:
        newat.append(np.dot(rotmat,at))
    return newat

# Inner function
def findZaxis_SOM_rotation(atoms, angle):
    """
    Find axis minimizing CSM under n-fold rotation.
    
    Uses random search (900 trials) + Powell optimization to locate
    rotation axis where CSM(atoms, R(angle)·atoms) is minimal.
    
    Args:
        atoms: Ligand positions (Nx3)
        angle: Rotation angle (π/2 for 4-fold, 2π/3 for 3-fold, etc.)
    
    Returns:
        tuple: (axis_vector, csm_value)
            axis_vector: Best-fit rotation axis
            csm_value: Achieved CSM (< 1 indicates good symmetry)
    """
    def fitfun(xx):
        theta,phi = xx
        return ContinuousShapeMeasure(atoms, rotateArbAxis(atoms, theta, phi, angle))

    ## Select starting parameters
    startX0s = []
    startFF = []
    for i in range(900):
        x,y,z = np.random.uniform(-1,1,size=3)
        norm = np.sqrt((x**2 + y**2 + z**2))
        if norm <= 1:
            startX0s.append([np.arcsin(z/norm), np.arctan(y/x)])
            startFF.append(fitfun(startX0s[-1]))        
    x0 = startX0s[np.argmin(startFF)]

    res = minimize(fitfun, x0=x0, method='Powell')
    return np.array(anglesToVector(*res.x)), res.fun



