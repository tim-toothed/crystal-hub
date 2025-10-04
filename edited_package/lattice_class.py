"""
Crystallographic lattice utilities for unit cell manipulations.

Handles conversions between different coordinate systems in crystallography:
- Real space (Cartesian XYZ)
- Fractional coordinates (ABC crystal axes)
- Reciprocal space (inverse Angstroms)

Provides unit cell construction from lattice parameters and coordinate transformations
for crystal structure calculations.
"""

import numpy as np

class lattice:
	"""
    Crystallographic unit cell with coordinate transformation methods.
    
    Constructs real and reciprocal lattice vectors from lattice parameters
    (a, b, c, α, β, γ) and provides conversions between coordinate systems.
    
    Attributes:
        a, b, c: Real space lattice vectors in Cartesian coordinates (Å)
        anorm, bnorm, cnorm: Lattice parameter lengths (Å)
        astar, bstar, cstar: Reciprocal lattice vectors (2π/Å)
        V: Unit cell volume (Å³)
        Vstar: Reciprocal cell volume
    
    Coordinate systems:
        - Cartesian: Standard XYZ coordinates
        - ABC: Fractional coordinates along crystal axes
        - RLU: Reciprocal lattice units (for diffraction)
    """

	def __init__(self, length_a, length_b, length_c, angle_alpha, angle_beta, angle_gamma):
		"""
        Initialize unit cell from lattice parameters.
        
        Args:
            length_a, length_b, length_c: Lattice parameters in Angstroms
            angle_alpha, angle_beta, angle_gamma: Lattice angles in degrees
                alpha: angle between b and c
                beta: angle between a and c
                gamma: angle between a and b
		"""

		self.anorm = length_a
		self.bnorm = length_b
		self.cnorm = length_c
		angle_alpha = angle_alpha*np.pi/180
		angle_beta = angle_beta*np.pi/180
		angle_gamma = angle_gamma*np.pi/180
		self.a = length_a*np.array([1,0,0])
		self.b = length_b*np.array([np.cos(angle_gamma), np.sin(angle_gamma), 0])
		self.c = length_c*np.array([np.cos(angle_beta), 
			np.sin(angle_beta)*np.cos(angle_alpha)*np.sin(angle_gamma), 
			np.sin(angle_alpha)*np.sin(angle_beta) ])
		#Round numbers (eliminate tiny numbers)
		self.b=np.around(self.b,decimals=8)
		self.c=np.around(self.c,decimals=8)

		#Define the reciprocal lattice
		self.reciplatt()

	def reciplatt(self):
		"""
        Calculate reciprocal lattice vectors.
        
        Computes a*, b*, c* from real space vectors using:
            a* = 2π (b × c) / V
        where V is the unit cell volume.
        
        Sets attributes: astar, bstar, cstar, V, Vstar
        """

		cellvol = np.dot(self.a,np.cross(self.b,self.c))

		self.astar = 2*np.pi * np.cross(self.b,self.c) / cellvol
		self.bstar = 2*np.pi * np.cross(self.c,self.a) / cellvol
		self.cstar = 2*np.pi * np.cross(self.a,self.b) / cellvol

		#Round numbers (eliminate tiny numbers)
		self.bstar=np.around(self.bstar,decimals=8)
		self.cstar=np.around(self.cstar,decimals=8)

		self.V = cellvol
		self.Vstar = np.dot(self.astar,np.cross(self.bstar,self.cstar))

	def cartesian(self,vect,norm=True):
		"""
        Convert from fractional (ABC) to Cartesian (XYZ) coordinates.
        
        Args:
            vect: 3-component vector or array of vectors in ABC coordinates
            norm: If True, use actual lattice vectors. If False, use normalized vectors.
        
        Returns:
            Vector(s) in Cartesian coordinates (Angstroms)
        
        Example:
            >>> latt = lattice(5, 5, 5, 90, 90, 90)  # Cubic
            >>> cart = latt.cartesian([0.5, 0.5, 0.5])  # Convert [½,½,½] to XYZ
		"""

		if vect.size == 3:
			if norm == True:
				return vect[0]*self.a + vect[1]*self.b + vect[2]*self.c
			else:
				return vect[0]*self.a/np.linalg.norm(self.a) \
					+ vect[1]*self.b/np.linalg.norm(self.b) \
					+ vect[2]*self.c/np.linalg.norm(self.c)
		elif len(vect[0]) == 3:
			return np.outer(vect[:,0],self.a) + np.outer(vect[:,1],self.b) + np.outer(vect[:,2],self.c)
		else:
			raise ValueError("vector must have three components") 

	def ABC(self, vect):
		"""
        Convert from Cartesian (XYZ) to fractional (ABC) coordinates.
        
        Args:
            vect: 3-component vector in Cartesian coordinates (Angstroms)
        
        Returns:
            Vector in fractional ABC coordinates
        
        Raises:
            ValueError: If vector doesn't have 3 components
        """

		matrix = np.array([self.a, self.b, self.c])

		if vect.size == 3:
			return np.dot(vect,np.linalg.inv(matrix))
		else:
			raise ValueError("vector must have three components") 

	def inverseA(self,vect,norm=True):
		"""
        Convert from reciprocal lattice units (RLU) to inverse Angstroms.
        
        Transforms momentum transfer vectors from Miller indices to absolute
        reciprocal space coordinates.
        
        Args:
            vect: 3-component vector or array in RLU (h,k,l)
            norm: Currently unused parameter (kept for compatibility)
        
        Returns:
            Vector(s) in inverse Angstroms
        
        Example:
            >>> Q = latt.inverseA([1, 0, 0])  # Convert (100) reflection
		"""
		if vect.size == 3:
			return vect[0]*self.astar + vect[1]*self.bstar + vect[2]*self.cstar
		elif len(vect[0])==3:
			return np.outer(vect[:,0],self.astar) + np.outer(vect[:,1],self.bstar) + np.outer(vect[:,2],self.cstar)
		else:
			raise ValueError("vector must have three components") 
