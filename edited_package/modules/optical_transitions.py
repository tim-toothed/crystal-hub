"""
Numba-optimized class for calculating optical transition probabilities.
"""

from numba import jitclass, float64
import numpy as np

# Specification for Numba jitclass - defines the data types
_spec = [ 
    ('Jx', float64[:,:]),  # J_x operator matrix
    ('Jy', float64[:,:]),  # J_y operator matrix
    ('Jz', float64[:,:])   # J_z operator matrix
]

@jitclass(_spec)
class OpticalTransition:
    """
    Numba-optimized class for calculating transition probabilities between 
    quantum states using angular momentum operators.
    
    Attributes:
        Jx (ndarray): x-component angular momentum operator matrix
        Jy (ndarray): y-component angular momentum operator matrix  
        Jz (ndarray): z-component angular momentum operator matrix
    """
    
    def __init__(self, Jx_matrix, Jy_matrix, Jz_matrix):
        """
        Initialize with angular momentum operator matrices.
        
        Args:
            Jx_matrix (ndarray): Matrix representation of J_x operator
            Jy_matrix (ndarray): Matrix representation of J_y operator
            Jz_matrix (ndarray): Matrix representation of J_z operator
        """
        self.Jx = Jx_matrix
        self.Jy = Jy_matrix
        self.Jz = Jz_matrix

    def transition_strength(self, bra_state, ket_state):
        """
        Calculate the total transition strength between two quantum states.
        
        The transition strength is proportional to |<bra|J|ket>|², summed over
        all three spatial components (x, y, z).
        
        Args:
            bra_state (ndarray): Bra vector (complex conjugate of final state)
            ket_state (ndarray): Ket vector (initial state)
            
        Returns:
            float: Total transition strength between the two states
        """
        # Calculate matrix elements for each component
        Jx_element = np.dot(bra_state, np.dot(self.Jx, ket_state))
        Jy_element = np.dot(bra_state, np.dot(self.Jy, ket_state)) 
        Jz_element = np.dot(bra_state, np.dot(self.Jz, ket_state))
        
        # Sum squared magnitudes for total transition strength
        return (np.abs(Jx_element)**2 + 
                np.abs(Jy_element)**2 + 
                np.abs(Jz_element)**2)

    # Alias for backward compatibility
    transition = transition_strength
