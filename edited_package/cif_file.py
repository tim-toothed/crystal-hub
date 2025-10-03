"""
Simple CIF (Crystallographic Information File) Importer
For crystal structure analysis in physics applications.

This module provides functionality to:
- Parse CIF files and extract crystallographic data
- Generate unit cells from asymmetric units using symmetry operations
- Calculate nuclear structure factors
- Analyze multiple scattering paths
"""

import numpy as np
from copy import deepcopy
from .lattice_class import Lattice


class CifFile:
    """Simple and efficient CIF file importer for crystal structure analysis."""
    
    def __init__(self, filepath):
        """Initialize and import CIF file."""
        self.filepath = filepath
        self.asymunitcell = []
        self.unitcell = []
        self.symops = []
        self.atomsunitcell = {}
        
        # Parse the file
        self._parse_cif()
        
        # Create Lattice object
        self.latt = Lattice(
            self.a, self.b, self.c,
            self.alpha, self.beta, self.gamma
        )
        
        print("CIF import complete.")
    
    def _parse_cif(self):
        """Parse CIF file and extract all necessary data."""
        with open(self.filepath, 'r') as f:
            lines = f.readlines()
        
        i = 0
        sites = []
        symops = []
        data_blocks = 0
        
        while i < len(lines):
            line = lines[i]
            
            # Check for phase boundaries
            if 'data_global' in line:
                data_blocks += 1
                if data_blocks > 1:
                    break
            if '=END' in line:
                break
            
            # Parse unit cell parameters
            if line.startswith('_cell_length_a'):
                self.a = self._parse_number(line.split()[1])
            elif line.startswith('_cell_length_b'):
                self.b = self._parse_number(line.split()[1])
            elif line.startswith('_cell_length_c'):
                self.c = self._parse_number(line.split()[1])
            elif line.startswith('_cell_angle_alpha'):
                self.alpha = self._parse_number(line.split()[1])
            elif line.startswith('_cell_angle_beta'):
                self.beta = self._parse_number(line.split()[1])
            elif line.startswith('_cell_angle_gamma'):
                self.gamma = self._parse_number(line.split()[1])
                print(f'Unit cell: {self.a}, {self.b}, {self.c}, {self.alpha}, {self.beta}, {self.gamma}')
            
            # Parse atomic positions
            elif line.startswith("loop_") and i + 1 < len(lines):
                next_line = lines[i + 1].strip()
                if next_line.startswith("_atom_site_label") or next_line.startswith("_atom_site_type"):
                    i = self._parse_atoms(lines, i + 1, sites)
            
            # Parse symmetry operations
            elif line.startswith("loop_") and i + 1 < len(lines):
                if any(keyword in lines[i + 1] for keyword in [
                    "_space_group_symop", "_symmetry_equiv_pos"
                ]):
                    i = self._parse_symmetry(lines, i + 1, symops)
            
            i += 1
        
        # Validate results
        if not sites:
            raise RuntimeError("No atomic sites were found in the CIF file")
        
        if not symops:
            print("Warning: No symmetry operations found, using P1")
            symops = ['x,y,z']
        
        self.asymunitcell = sites
        self.symops = symops
        
        # Generate unit cell
        self.MakeUnitCell(sites, symops)
    
    def _parse_atoms(self, lines, start_idx, sites):
        """Parse atomic positions from CIF loop."""
        print("Importing atoms")
        i = start_idx
        
        # First, identify field indices
        fields = {}
        field_count = 0
        
        while i < len(lines) and lines[i].strip().startswith('_atom'):
            line = lines[i].strip()
            if 'atom_site_fract_x' in line:
                fields['x'] = field_count
            elif 'atom_site_fract_y' in line:
                fields['y'] = field_count
            elif 'atom_site_fract_z' in line:
                fields['z'] = field_count
            elif 'atom_site_occupancy' in line:
                fields['occ'] = field_count
            elif 'atom_site_site_symmetry_order' in line or 'atom_site_symmetry_multiplicity' in line:
                fields['symorder'] = field_count
            elif 'site_type_symbol' in line:
                fields['element'] = field_count
            elif 'site_label' in line:
                fields['label'] = field_count
            field_count += 1
            i += 1
        
        # Parse atom data
        while i < len(lines):
            line = lines[i].strip()
            if not line or line.startswith('loop_') or line.startswith('#') or line.startswith('_'):
                break
            
            parts = line.split()
            if len(parts) < 3:
                i += 1
                continue
            
            # Extract atom data
            atom = [''] * 10  # Pre-allocate list
            atom[0] = parts[fields.get('label', 0)]
            atom[1] = parts[fields.get('element', 1)]
            
            # Get fractional coordinates and wrap to unit cell
            x = self._parse_number(parts[fields['x']])
            y = self._parse_number(parts[fields['y']])
            z = self._parse_number(parts[fields['z']])
            
            # Wrap to [0, 1)
            atom[2] = x - int(x)
            if atom[2] < 0: atom[2] += 1
            atom[3] = y - int(y)
            if atom[3] < 0: atom[3] += 1
            atom[4] = z - int(z)
            if atom[4] < 0: atom[4] += 1
            
            # Padding for compatibility
            atom[5] = atom[6] = atom[7] = 0
            
            # Occupancy
            atom[7] = self._parse_number(parts[fields['occ']]) if 'occ' in fields else 1.0
            
            # Symmetry order
            atom[8] = int(parts[fields['symorder']]) if 'symorder' in fields else 0
            
            # Store label at end for compatibility
            atom[9] = atom[0]
            
            sites.append(atom)
            i += 1
        
        return i - 1
    
    def _parse_symmetry(self, lines, start_idx, symops):
        """Parse symmetry operations from CIF loop."""
        i = start_idx
        
        # Skip header lines
        while i < len(lines) and '_' in lines[i]:
            i += 1
        
        # Parse operations
        while i < len(lines):
            line = lines[i].strip()
            if not line or line.startswith('loop_') or line.startswith('#'):
                break
            
            # Extract symmetry operation
            if '\'' in line:
                quotes = [j for j, ch in enumerate(line) if ch == '\'']
                if len(quotes) >= 2:
                    symops.append(line[quotes[0]+1:quotes[1]])
            else:
                for part in line.split():
                    if ',' in part:
                        symops.append(part.strip())
                        break
            i += 1
        
        return i - 1
    
    def _parse_number(self, string):
        """Parse numeric value from CIF, handling uncertainties and fractions."""
        if '(' in string:
            string = string[:string.index('(')]
        
        value = float(string)
        
        # Check for common fractions
        threshold = 5e-5
        if abs(value - 1/3) < threshold:
            return 1/3
        elif abs(value - 2/3) < threshold:
            return 2/3
        elif abs(value - 1/6) < threshold:
            return 1/6
        elif abs(value - 5/6) < threshold:
            return 5/6
        
        return value
    
    def apply_sym_operation(self, symstring, atom):
        """Apply symmetry operation to atomic position."""
        newatom = list(atom)
        x, y, z = atom[2], atom[3], atom[4]
        
        # Clean and substitute coordinates
        symstring = symstring.replace(' ', '')
        symstring = (symstring
                    .replace('x', str(x)).replace('X', str(x))
                    .replace('y', str(y)).replace('Y', str(y))
                    .replace('z', str(z)).replace('Z', str(z)))
        
        # Evaluate each component
        components = symstring.split(',')
        newatom[2] = eval(components[0])
        newatom[3] = eval(components[1])
        newatom[4] = eval(components[2])
        
        # Wrap back to unit cell
        for i in range(2, 5):
            if newatom[i] < -0.001:
                newatom[i] += 1
            if newatom[i] >= 1.0:
                newatom[i] -= 1
        
        return newatom
    
    def MakeUnitCell(self, sites, symops):
        """Generate complete unit cell from asymmetric unit and symmetry operations."""
        unitcell = []
        
        # Apply all symmetry operations to all sites
        for symop in symops:
            for site in sites:
                new_atom = self.apply_sym_operation(symop, site)
                
                # Wrap to unit cell
                for i in range(2, 5):
                    new_atom[i] = new_atom[i] - int(new_atom[i])
                    if new_atom[i] < 0:
                        new_atom[i] += 1
                
                # Check for duplicates
                is_duplicate = False
                for existing in unitcell:
                    if (existing[1] == new_atom[1] and  # Same element
                        abs(existing[2] - new_atom[2]) < 0.001 and
                        abs(existing[3] - new_atom[3]) < 0.001 and
                        abs(existing[4] - new_atom[4]) < 0.001):
                        is_duplicate = True
                        break
                
                if not is_duplicate:
                    unitcell.append(new_atom)
        
        # Keep only atoms inside unit cell
        self.unitcell = []
        for atom in unitcell:
            if (0 <= atom[2] <= 1 and 
                0 <= atom[3] <= 1 and 
                0 <= atom[4] <= 1):
                self.unitcell.append(atom)
        
        print(f'  {len(self.unitcell)} atoms added')
        
        # Create atom position dictionary by label
        self.atomsunitcell = {}
        for atom in self.unitcell:
            label = atom[9]  # Label stored at end
            if label not in self.atomsunitcell:
                self.atomsunitcell[label] = []
            self.atomsunitcell[label].append([atom[2], atom[3], atom[4]])
        
        # Convert to numpy arrays
        for label in self.atomsunitcell:
            self.atomsunitcell[label] = np.array(self.atomsunitcell[label])
    
    def StructureFactor(self, scattering_lengths, max_hkl):
        """
        Calculate nuclear structure factors.
        
        Args:
            scattering_lengths: Dict of neutron scattering lengths by element (fm)
            max_hkl: Maximum Miller index to calculate
        
        Returns:
            Array of [h, k, l, |F|^2] for each reflection
        """
        # Validate scattering lengths
        elements = set(atom[1] for atom in self.unitcell)
        missing = [el for el in elements if el not in scattering_lengths]
        if missing:
            raise ValueError(
                f"Missing scattering lengths for: {missing}\n"
                "See: https://www.ncnr.nist.gov/resources/n-lengths/"
            )
        
        # Build atom array: [x, y, z, b_coherent, occupancy]
        atom_data = np.zeros((len(self.unitcell), 5))
        for i, atom in enumerate(self.unitcell):
            atom_data[i] = [
                atom[2], atom[3], atom[4],
                scattering_lengths[atom[1]],
                atom[7]  # occupancy
            ]
        
        # Generate reflection list
        hkl_range = np.arange(-max_hkl + 1, max_hkl)
        h, k, l = np.meshgrid(hkl_range, hkl_range, hkl_range)
        reflections = np.column_stack([h.ravel(), k.ravel(), l.ravel()])
        
        # Calculate structure factors
        sf = np.zeros((len(reflections), 4))
        sf[:, :3] = reflections
        
        for i, hkl in enumerate(reflections):
            # F = Σ b_j * occ_j * exp(2πi(h·r_j))
            phase = 2j * np.pi * np.dot(atom_data[:, :3], hkl)
            F = np.sum(atom_data[:, 3] * atom_data[:, 4] * np.exp(phase))
            sf[i, 3] = abs(F) ** 2
        
        # Remove tiny values
        sf[:, 3][sf[:, 3] < 1e-16] = 0.0
        
        self.SF = sf
        return sf
    
    def MultipleScattering(self, ei, peak, xcut, ycut, threshold=0.05):
        """
        Analyze multiple scattering paths for a given reflection.
        
        Args:
            ei: Incident neutron energy (meV)
            peak: Target reflection [h, k, l]
            xcut: First vector defining scattering plane
            ycut: Second vector defining scattering plane  
            threshold: Maximum deviation in Q-space
        """
        if not hasattr(self, 'SF'):
            raise RuntimeError("Must calculate structure factors first")
        
        # Calculate incident wavevector (Å^-1)
        k = 0.694693 * np.sqrt(ei)
        
        # Normalize scattering plane vectors
        xcut_norm = xcut / np.linalg.norm(xcut)
        ycut_norm = ycut / np.linalg.norm(ycut)
        
        # Calculate scattering angle for target peak
        q_peak = np.linalg.norm(self.latt.inverseA(vect=np.array(peak)))
        theta = np.arcsin(min(q_peak / (2 * k), 1.0))
        
        # Incident k in reciprocal space
        k_h = k * (xcut_norm[0] * np.sin(theta) + ycut_norm[0] * np.cos(theta))
        k_k = k * (xcut_norm[1] * np.sin(theta) + ycut_norm[1] * np.cos(theta))
        k_l = k * (xcut_norm[2] * np.sin(theta) + ycut_norm[2] * np.cos(theta))
        
        print(f"\nMultiple scattering for {peak} peak, Ei = {ei} meV:")
        print("-" * 45)
        print("    dq \t \t intermediate peaks \t \t  SF^2")
        
        # Find multiple scattering paths
        for t1 in self.SF:
            if t1[3] == 0 or np.allclose(t1[:3], [0, 0, 0]):
                continue
            
            # Calculate momentum transfer
            t1_A = self.latt.inverseA(vect=t1[:3])
            dq = np.sqrt((t1_A[0] - k_h)**2 + (t1_A[1] - k_k)**2 + (t1_A[2] - k_l)**2) - k
            
            if abs(dq) < threshold:
                # Find second scattering
                t2_hkl = np.array(peak) - t1[:3]
                
                # Get structure factor for second reflection
                t2_sf = 0
                for refl in self.SF:
                    if np.allclose(refl[:3], t2_hkl):
                        t2_sf = refl[3]
                        break
                
                if t2_sf > 0 and not np.allclose(t2_hkl, [0, 0, 0]):
                    # Check scattering angles are reasonable
                    q1 = np.linalg.norm(t1[:3]) * 2 * np.pi / 10.0
                    q2 = np.linalg.norm(t2_hkl) * 2 * np.pi / 10.0
                    
                    if q1 < 2*k and q2 < 2*k:
                        total_sf = t1[3] * t2_sf
                        print(f"   {dq:6.4f}\t {t1[:3].astype(int)}  {t2_hkl.astype(int)} \t {total_sf}")


# Example usage
if __name__ == "__main__":
    # Test with Yb2Ti2O7
    crystal = CifFile("StoichiometricYbTiO.cif")
    
    # Define scattering lengths
    s_length = {'O2-': 5.803, 'Ti4+': -3.438, 'Yb3+': 12.43}
    
    # Calculate structure factors
    crystal.StructureFactor(s_length, 5)
    print(f"\nCalculated {len(crystal.SF)} reflections")
    
    # Analyze multiple scattering
    crystal.MultipleScattering(
        ei=1.0, 
        peak=[0, 0, 2],
        xcut=np.array([1, 1, 1]),
        ycut=np.array([1, 1, -2]),
        threshold=0.1
    )