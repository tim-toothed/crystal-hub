"""
cif_parser.py - CIF/mCIF to NJA-CFS converter

Converts CIF and mCIF (magnetic CIF) files to the tab-delimited format required by NJA-CFS
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Union
import warnings

def parse_cif_to_nja(
    cif_file: str,
    metal_center: Union[str, int],
    coordination_cutoff: float = 3.5,
    charge_dict: Optional[Dict[str, float]] = None,
    output_file: Optional[str] = None,
    use_symmetry: bool = True,
    radial_displacement: float = 0.0,
    file_type: str = 'auto'
) -> np.ndarray:
    """
    Parse CIF or mCIF file and convert to NJA-CFS input format.
    
    Parameters:
    -----------
    cif_file : str
        Path to CIF or mCIF file
    metal_center : str or int
        Element symbol (e.g., 'Tb') or atom index of metal center
    coordination_cutoff : float
        Maximum distance (Å) to include ligands
    charge_dict : dict, optional
        Dictionary mapping element symbols to charges
        Example: {'O': -2.0, 'N': -1.0, 'C': 0.0}
        If None, uses default oxidation states
    output_file : str, optional
        If provided, saves to this file
    use_symmetry : bool
        Whether to apply symmetry operations to generate full coordination sphere
    radial_displacement : float
        Radial displacement parameter (usually 0.0)
    file_type : str
        'cif', 'mcif', or 'auto' to auto-detect
    
    Returns:
    --------
    data : np.ndarray
        Array with columns: [label, x, y, z, charge]
    """
    
    try:
        from pymatgen.io.cif import CifParser
        from pymatgen.core import Structure
    except ImportError:
        raise ImportError(
            "pymatgen is required for CIF parsing. Install with:\n"
            "pip install pymatgen"
        )
    
    # Auto-detect file type
    if file_type == 'auto':
        if cif_file.lower().endswith('.mcif'):
            file_type = 'mcif'
            print("Detected mCIF (magnetic CIF) file")
        else:
            file_type = 'cif'
            print("Detected CIF file")
    
    # Default charge assignments
    default_charges = {
        'O': -2.0, 'S': -2.0, 'Se': -2.0, 'Te': -2.0,
        'N': -1.0, 'P': -1.0, 'As': -1.0,
        'F': -1.0, 'Cl': -1.0, 'Br': -1.0, 'I': -1.0,
        'C': 0.0, 'H': 0.0, 'B': 0.0,
        # Add common metals with their typical oxidation states
        'Ti': 4.0, 'V': 3.0, 'Cr': 3.0, 'Mn': 2.0, 
        'Fe': 3.0, 'Co': 2.0, 'Ni': 2.0, 'Cu': 2.0, 'Zn': 2.0,
        # Lanthanides (typically 3+)
        'La': 3.0, 'Ce': 3.0, 'Pr': 3.0, 'Nd': 3.0, 'Pm': 3.0,
        'Sm': 3.0, 'Eu': 3.0, 'Gd': 3.0, 'Tb': 3.0, 'Dy': 3.0,
        'Ho': 3.0, 'Er': 3.0, 'Tm': 3.0, 'Yb': 3.0, 'Lu': 3.0,
    }
    
    if charge_dict is None:
        charge_dict = default_charges
        warnings.warn(
            "No charge_dict provided. Using default oxidation states. "
            "For accurate results, provide explicit charges for your system.",
            UserWarning
        )
    
    # Parse file (both CIF and mCIF are handled the same way by pymatgen)
    try:
        parser = CifParser(cif_file)
        # Try new API first, fall back to old if needed
        try:
            structure = parser.parse_structures(primitive=False)[0]
        except AttributeError:
            structure = parser.get_structures()[0]
    except Exception as e:
        raise ValueError(f"Error parsing file: {e}\n"
                        f"Make sure the file is a valid CIF or mCIF format.")
    
    print(f"\n{'='*70}")
    print(f"Structure Information")
    print(f"{'='*70}")
    print(f"Formula: {structure.composition.reduced_formula}")
    try:
        print(f"Space group: {structure.get_space_group_info()}")
    except:
        print("Space group: Not available")
    
    print(f"\nUnit cell parameters:")
    print(f"  a = {structure.lattice.a:.4f} Å")
    print(f"  b = {structure.lattice.b:.4f} Å")
    print(f"  c = {structure.lattice.c:.4f} Å")
    print(f"  α = {structure.lattice.alpha:.2f}°")
    print(f"  β = {structure.lattice.beta:.2f}°")
    print(f"  γ = {structure.lattice.gamma:.2f}°")
    
    print(f"\nAtoms in asymmetric unit:")
    for i, site in enumerate(structure):
        print(f"  {i}: {site.species_string:8s} at fractional coords {site.frac_coords}")
    
    # Find metal center
    if isinstance(metal_center, str):
        from pymatgen.core import Element
        
        # Try to get element from string
        try:
            target_element = Element(metal_center.rstrip('0123456789+-'))
        except:
            target_element = Element(metal_center)
        
        metal_indices = []
        for i, site in enumerate(structure):
            site_elements = [el.symbol for el in site.species.elements]
            if target_element.symbol in site_elements:
                metal_indices.append(i)
        
        if not metal_indices:
            available_elements = set(el.symbol for site in structure 
                                    for el in site.species.elements)
            raise ValueError(
                f"No {metal_center} atoms found in structure.\n"
                f"Available elements: {sorted(available_elements)}"
            )
        
        if len(metal_indices) > 1:
            print(f"\n{'='*70}")
            print(f"Found {len(metal_indices)} {target_element.symbol} sites:")
            print(f"{'='*70}")
            for idx in metal_indices:
                site = structure[idx]
                occupancy = site.species.get_atomic_fraction(target_element)
                print(f"  Index {idx}: {site.species_string:8s} "
                      f"(occupancy: {occupancy:.3f}) at {site.frac_coords}")
            
            # Use the one with highest occupancy
            occupancies = [structure[idx].species.get_atomic_fraction(target_element) 
                          for idx in metal_indices]
            metal_idx = metal_indices[np.argmax(occupancies)]
            print(f"\nUsing site with highest occupancy: index {metal_idx}")
        else:
            metal_idx = metal_indices[0]
    else:
        metal_idx = int(metal_center)
        if metal_idx >= len(structure):
            raise ValueError(f"Index {metal_idx} out of range. "
                           f"Structure has {len(structure)} sites.")
    
    metal_site = structure[metal_idx]
    metal_coords = metal_site.coords
    
    print(f"\n{'='*70}")
    print(f"Selected metal center:")
    print(f"{'='*70}")
    print(f"  Element: {metal_site.species_string}")
    print(f"  Cartesian coords: {metal_coords}")
    print(f"  Fractional coords: {metal_site.frac_coords}")
    
    # Find coordination sphere
    neighbors = structure.get_neighbors(metal_site, coordination_cutoff)
    
    if not neighbors:
        raise ValueError(
            f"No neighbors found within {coordination_cutoff} Å.\n"
            f"Try increasing coordination_cutoff parameter."
        )
    
    print(f"\nFound {len(neighbors)} atoms within {coordination_cutoff} Å")
    
    # Build data array
    data_list = []
    
    for neighbor in neighbors:
        site = neighbor[0]  # neighbor is tuple (site, distance, index)
        distance = neighbor[1]
        
        # Get element symbol (strip oxidation state and numbers)
        element = str(site.specie).rstrip('0123456789+-')
        
        # Get charge
        if element not in charge_dict:
            warnings.warn(
                f"Element {element} not in charge_dict. Assigning charge 0.0.\n"
                f"Consider adding '{element}': <charge> to charge_dict for accurate results.",
                UserWarning
            )
            charge = 0.0
        else:
            charge = charge_dict[element]
        
        # Calculate relative coordinates (metal at origin)
        rel_coords = site.coords - metal_coords
        
        # Create unique label (element + sequential number)
        label = f"{element}{len([d for d in data_list if d[0].startswith(element)]) + 1}"
        
        data_list.append([
            label,
            rel_coords[0],
            rel_coords[1],
            rel_coords[2],
            charge
        ])
    
    # Convert to numpy array
    data = np.array(data_list, dtype=object)
    
    # Sort by distance from metal center
    distances = np.sqrt(
        data[:, 1].astype(float)**2 + 
        data[:, 2].astype(float)**2 + 
        data[:, 3].astype(float)**2
    )
    sorted_indices = np.argsort(distances)
    data = data[sorted_indices]
    
    # Print coordination sphere summary
    print(f"\n{'='*70}")
    print(f"Coordination Sphere (sorted by distance)")
    print(f"{'='*70}")
    print(f"{'Label':<8} {'x (Å)':>9} {'y (Å)':>9} {'z (Å)':>9} {'d (Å)':>9} {'Charge':>8}")
    print("-" * 70)
    for i, row in enumerate(data):
        d = distances[sorted_indices[i]]
        print(f"{row[0]:<8} {float(row[1]):9.4f} {float(row[2]):9.4f} "
              f"{float(row[3]):9.4f} {d:9.4f} {float(row[4]):8.2f}")
    
    # Calculate coordination number
    print(f"\nCoordination number: {len(data)}")
    
    # Save to file if requested
    if output_file:
        save_nja_input(data, output_file, radial_displacement)
        print(f"\nSaved to: {output_file}")
    
    return data


def save_nja_input(
    data: np.ndarray,
    filename: str,
    radial_displacement: float = 0.0
):
    """
    Save data in NJA-CFS tab-delimited format.
    
    Parameters:
    -----------
    data : np.ndarray
        Array with columns: [label, x, y, z, charge]
    filename : str
        Output filename
    radial_displacement : float
        Radial displacement parameter
    """
    with open(filename, 'w') as f:
        f.write("# NJA-CFS input file\n")
        f.write("# Format: label  charge  r_displacement  x  y  z\n")
        f.write("# Coordinates in Angstroms, metal center at origin\n")
        f.write("#\n")
        for row in data:
            label = row[0]
            charge = float(row[4])
            x = float(row[1])
            y = float(row[2])
            z = float(row[3])
            
            f.write(f"{label}\t{charge:.6f}\t{radial_displacement:.6f}\t"
                   f"{x:.6f}\t{y:.6f}\t{z:.6f}\n")


def cif_to_bkq(
    cif_file: str,
    metal_center: Union[str, int],
    configuration: str,
    coordination_cutoff: float = 3.5,
    charge_dict: Optional[Dict[str, float]] = None,
    use_sternheimer: bool = False,
    file_type: str = 'auto'
) -> Dict:
    """
    Direct conversion from CIF/mCIF to Bkq parameters.
    
    Parameters:
    -----------
    cif_file : str
        Path to CIF or mCIF file
    metal_center : str or int
        Element symbol or atom index of metal center
    configuration : str
        Electron configuration (e.g., 'f7', 'd5')
    coordination_cutoff : float
        Maximum distance (Å) to include ligands
    charge_dict : dict, optional
        Dictionary mapping element symbols to charges
    use_sternheimer : bool
        Whether to apply Sternheimer shielding parameters
    file_type : str
        'cif', 'mcif', or 'auto' to auto-detect
    
    Returns:
    --------
    dic_bkq : dict
        Crystal field parameters in Wybourne notation
    """
    
    import nja_cfs_red as nja
    
    # Parse CIF/mCIF
    data = parse_cif_to_nja(
        cif_file,
        metal_center,
        coordination_cutoff,
        charge_dict,
        file_type=file_type
    )
    
    # Calculate Bkq parameters
    dic_bkq = nja.calc_Bkq(
        data,
        configuration,
        sph_flag=False,
        sth_param=use_sternheimer
    )
    
    return dic_bkq


def visualize_coordination_sphere(
    data: np.ndarray,
    metal_label: str = "Metal",
    save_fig: Optional[str] = None,
    show_charges: bool = True
):
    """
    Create 3D visualization of coordination sphere.
    
    Parameters:
    -----------
    data : np.ndarray
        Array with columns: [label, x, y, z, charge]
    metal_label : str
        Label for metal center
    save_fig : str, optional
        If provided, saves figure to this file
    show_charges : bool
        Whether to color atoms by charge
    """
    try:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
    except ImportError:
        print("matplotlib required for visualization")
        return
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot metal center
    ax.scatter(0, 0, 0, c='gold', s=300, marker='o', 
              edgecolors='black', linewidths=2, label=metal_label, zorder=10)
    
    # Plot ligands
    coords = data[:, 1:4].astype(float)
    charges = data[:, 4].astype(float)
    labels = data[:, 0]
    
    if show_charges and len(np.unique(charges)) > 1:
        # Color by charge
        norm = plt.cm.colors.Normalize(vmin=charges.min(), vmax=charges.max())
        colors = plt.cm.RdBu(norm(charges))
        scatter = ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
                           c=charges, cmap='RdBu', s=150, marker='o', 
                           alpha=0.8, edgecolors='black', linewidths=1)
        plt.colorbar(scatter, ax=ax, label='Charge', shrink=0.6)
    else:
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
                  c='lightblue', s=150, marker='o', alpha=0.8,
                  edgecolors='black', linewidths=1)
    
    # Add labels
    for i, (label, coord) in enumerate(zip(labels, coords)):
        ax.text(coord[0], coord[1], coord[2], f'  {label}', 
               fontsize=9, weight='bold')
    
    # Draw bonds
    for coord in coords:
        ax.plot([0, coord[0]], [0, coord[1]], [0, coord[2]], 
               'k-', alpha=0.2, linewidth=1)
    
    ax.set_xlabel('x (Å)', fontsize=11)
    ax.set_ylabel('y (Å)', fontsize=11)
    ax.set_zlabel('z (Å)', fontsize=11)
    ax.set_title(f'Coordination Sphere around {metal_label}', fontsize=13, weight='bold')
    
    # Equal aspect ratio
    max_range = np.abs(coords).max() * 1.1
    ax.set_xlim([-max_range, max_range])
    ax.set_ylim([-max_range, max_range])
    ax.set_zlim([-max_range, max_range])
    
    # Improve viewing angle
    ax.view_init(elev=20, azim=45)
    
    plt.tight_layout()
    
    if save_fig:
        plt.savefig(save_fig, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {save_fig}")
    
    plt.show()


if __name__ == "__main__":
    print(__doc__)