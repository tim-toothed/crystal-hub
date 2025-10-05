'''
<- подумать над тем, что это за модуль вообще
'''

import numpy as np
import matplotlib.pyplot as plt

# used outside the package

def plotPCF(onesite, nearestNeighbors, Xax, Yax, Zax, save_file: bool = False):
    """
    Plot Point Charge Field (PCF) ligand environment around a central ion.
    
    Creates a 2D projection showing the central ion at origin with surrounding
    ligand positions and connecting bonds. Coordinate axes are displayed.
    
    Args:
        onesite: Central ion information, onesite[0] is the ion name
        nearestNeighbors: List of nearest neighbor data, each entry has 
                         position at index [2]
        Xax, Yax, Zax: Coordinate system axis vectors to display
        save_file: If True, saves plot as '{ion_name}_ligands.png'. 
                  If False, displays interactively. Default False.
    
    Example:
        >>> plotPCF(central_ion, ligands, X, Y, Z, save_file=True)
    """
    atomlist = [[[0,0,0]], [nn[2] for nn in nearestNeighbors]]

    plt.figure(figsize=(6,6))
    obj = atomplot(1.5, 0.2, atomlist)
    obj.plotatoms(plotlines=[1])
    obj.plotaxes(Xax, Yax, Zax)
    
    if save_file:
        plt.savefig(onesite[0] + '_ligands.png')
        plt.close()
    else:
        plt.show()

# used inside the package
# <- Используется только для plotPCF
class atomplot:

	"""
    3D visualization tool for plotting atomic positions in crystal structures.
    
    Projects 3D atomic coordinates onto a 2D plane defined by viewing angles 
    (theta, phi). Used to visualize ligand environments around central ions.
    
    Attributes:
        theta: Polar angle for viewing direction (radians)
        phi: Azimuthal angle for viewing direction (radians)
        plotX, plotY, plotZ: Orthonormal basis vectors defining the viewing plane
        atoms: List of atomic position arrays to plot
    
    Methods:
        plotatoms(): Renders atoms as colored circles with optional bonds
        plotaxes(): Draws labeled X, Y, Z coordinate axes
        plotabc(): Draws crystallographic a, b, c axes in RGB colors
    """

	def __init__(self, theta, phi, atoms):
		"""
        Initialize the 3D atom plotter with viewing angles.
        
        Args:
            theta: Polar angle for viewing direction (radians)
            phi: Azimuthal angle for viewing direction (radians)
            atoms: List of atomic position groups. Each group is a list of 3D coordinates.
                   Example: [[[0,0,0]], [[1,0,0], [0,1,0]]] for central ion + 2 ligands
        """

		self.theta = theta
		self.phi = phi
		self.plotX = np.array([np.cos(phi), np.sin(phi), 0])
		self.plotY = np.array([-np.sin(phi)*np.cos(theta), np.cos(phi)*np.cos(theta), np.sin(theta)])
		self.plotZ = np.cross(self.plotX, self.plotY)

		self.atoms = atoms

	def plotatoms(self, plotlines  = []):
		"""
        Render atoms as colored circles with depth ordering and optional bonds.
        
        Projects 3D atom positions onto the viewing plane. Each atom group gets
        a different color. Atom size decreases with group index. Z-ordering
        creates depth perception.
        
        Args:
            plotlines: List of group indices to draw connecting lines between atoms.
                      Lines connect atoms roughly aligned with the central ion direction.
                      Example: [1] draws bonds for the second atom group (index 1)
        """

		plt.axis('off')
		colors = plt.cm.Set1(np.arange(8))
		for i,at in enumerate(self.atoms):
			for aa in at:
				plt.plot(np.dot(aa, self.plotX), np.dot(aa, self.plotY), 
					marker='o', markersize=100/(i+2), mec='k', color=colors[i], 
					zorder= np.dot(aa, self.plotZ))

			if i in plotlines:
				for a1 in at:
					dist = np.array(a1)-np.array(self.atoms[0][0])
					for a2 in at:
						vect = np.array(a1)-np.array(a2)
						if np.abs(np.dot(vect,  dist)) < np.dot(dist,dist)*1.2:
							plt.plot([np.dot(a1, self.plotX), np.dot(a2, self.plotX)], 
								[ np.dot(a1, self.plotY), np.dot(a2, self.plotY)], 
								color='grey', lw='3', zorder = np.mean([np.dot(a1, self.plotZ), np.dot(a2, self.plotZ)]))

		xvals = [np.dot(a, self.plotX) for a in at for at in self.atoms]
		yvals = [np.dot(a, self.plotY) for a in at for at in self.atoms]
		plt.xlim(np.min(xvals+yvals)*1.25, np.max(xvals+yvals)*1.25)
		plt.ylim(np.min(xvals+yvals)*1.25, np.max(xvals+yvals)*1.25)

	def plotaxes(self, X, Y, Z):
		"""
        Draw labeled coordinate system axes as black arrows.
        
        Args:
            X, Y, Z: 3D axis vectors to display. These are projected onto
                    the viewing plane and drawn as arrows from origin.
        """

		arrowatributes = {'head_width':0.06, 'overhang':0.1, 'color':'k'}

		plt.arrow(0,0, *self._flatten(X/2), **arrowatributes)
		plt.arrow(0,0, *self._flatten(Y/2), **arrowatributes)
		plt.arrow(0,0, *self._flatten(Z/2), **arrowatributes)

		disp = np.array([0.04,0.04])
		plt.text(*self._flatten(X/2)+disp, 'X')
		plt.text(*self._flatten(Y/2)+disp, 'Y')
		plt.text(*self._flatten(Z/2)+disp, 'Z')

	def plotabc(self):
		"""
        Draw crystallographic a, b, c axes in red, green, blue colors.
        
        Standard Cartesian axes ([1,0,0], [0,1,0], [0,0,1]) are drawn
        as colored arrows representing crystallographic directions.
        """

		X = np.array([1,0,0])
		Y = np.array([0,1,0])
		Z = np.array([0,0,1])
		arrowatributes = {'head_width':0.08, 'overhang':0.1}

		plt.arrow(0,0, *self._flatten(X/2), color='r', **arrowatributes)
		plt.arrow(0,0, *self._flatten(Y/2), color='g', **arrowatributes)
		plt.arrow(0,0, *self._flatten(Z/2), color='b', **arrowatributes)

		disp = np.array([0.04,0.04])
		plt.text(*self._flatten(X/2)+disp, 'X')
		plt.text(*self._flatten(Y/2)+disp, 'Y')
		plt.text(*self._flatten(Z/2)+disp, 'Z')

	def _flatten(self, vect):
		"""
        Project a 3D vector onto the 2D viewing plane.
        
        Args:
            vect: 3D vector to project
            
        Returns:
            Tuple (x, y) of 2D coordinates in the viewing plane
        """

		return np.dot(vect, self.plotX), np.dot(vect, self.plotY)




