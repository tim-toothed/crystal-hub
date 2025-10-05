"""
CIF file import and automatic point charge model generation.

Reads crystallographic information files (CIF), identifies magnetic ions and their
ligand environments, and constructs point charge model Hamiltonians. Handles both
rare earth (J-basis) and transition metal (LS-basis) systems automatically.

Main workflow:
    CIF → symmetry analysis → ligand detection → point charge model → CFLevels

Supports:
- Automatic ion detection from CIF
- Multiple crystallographic sites (A/B positions)
- Inversion symmetry detection
- Coordination sphere selection by distance/count
"""

import numpy as np
from copy import deepcopy

from .cif_symmetry_import import FindPointGroupSymOps
from .cif_file import CifFile
from .constants import ION_NUMS_RARE_EARTH, SPIN_ORBIT_COUPLING_CONSTANTS, RADIAL_INTEGRALS_TRANS_METAL, ION_HALF_FILLED, ION_NOT_HALF_FILLED
from .ligands import Ligands, LS_Ligands


def importCIF(ciffile, mag_ion = None, Zaxis = None, Yaxis = None, LS_Coupling = None,
                crystalImage=False, NumIonNeighbors=1, ForceImaginary=False, 
                ionL = None, ionS = None, CoordinationNumber = None, MaxDistance=None,
                ):
    """
    Import CIF and generate point charge crystal field model.
    
    Automatically detects ligand geometry, applies symmetry operations, and
    constructs Hamiltonian. Handles rare earth and transition metal ions.
    
    Args:
        ciffile: Path to CIF file
        mag_ion: Central ion label (e.g., 'Yb1', 'Ni1'). Auto-detected if None.
        Zaxis: Preferred z-axis direction for orientation
        Yaxis: Preferred y-axis direction
        LS_Coupling: Spin-orbit coupling λ (meV). Auto-detected for TM if None.
        crystalImage: If True, plot ligand geometry
        NumIonNeighbors: Include nth nearest neighbor shell
        ForceImaginary: If True, include m<0 terms (no inversion symmetry)
        ionL: Orbital angular momentum (required for TM ions not in database)
        ionS: Spin quantum number (required for TM ions not in database)
        CoordinationNumber: Force specific coordination (overrides auto-detection)
        MaxDistance: Maximum ligand distance cutoff (Angstroms)
    
    Returns:
        Single site: tuple (Ligands, CFLevels)
        Multiple sites: list [[Lig1, CF1], [Lig2, CF2], ...]
    
    Raises:
        TypeError: If TM ion requires ionL/ionS but not provided
    
    Example:
        >>> # Rare earth: auto-detect everything
        >>> lig, cf = importCIF('Yb2Ti2O7.cif')
        >>> cf.diagonalize()
        
        >>> # Transition metal: specify L, S
        >>> lig, cf = importCIF('NiO.cif', mag_ion='Ni1', ionL=2, ionS=1,
        ...                     LS_Coupling=-315)
    
    Note:
        - Handles multiple crystallographic sites (A/B, prime notation)
        - Automatically chooses J-basis (RE) vs LS-basis (TM)
        - Inversion symmetry detection sets suppressminusm flag
    """

    cif = CifFile(ciffile)
    if mag_ion == None: #take the first rare earth in the cif file as the central ion
        for asuc in cif.asymunitcell:
            if asuc[1].strip('3+') in ['Sm','Pm','Nd','Ce','Dy','Ho','Tm','Pr','Er','Tb','Yb']:
                mag_ion = asuc[0]
                print('No mag_ion ion listed, assuming', mag_ion, 'is the central ion.')
                break


    ## Check for multiply defined atoms
    differentPositionsA = []
    differentPositionsB = []
    for ii, at in enumerate(cif.unitcell):
        if at[4] < 0: print('negative atom!',ii, at)

        if at[0][-1] in ["'", "B", "b"]:
            differentPositionsA.append(at[0])
            differentPositionsB.append(at[0].replace("'","").replace("B","A").replace("b","a"))

    if len(differentPositionsA) > 0:
        cif_a = deepcopy(cif)
        cif_b = deepcopy(cif)

        unitcellA = []
        unitcellB = []
        for ii, at in enumerate(cif.unitcell):
            if at[0] in differentPositionsA:
                unitcellA.append(at)
            elif at[0] in differentPositionsB:
                unitcellB.append(at)
            else:
                unitcellA.append(at)
                unitcellB.append(at)

        cif_a.unitcell = unitcellA
        cif_b.unitcell = unitcellB

        cifs = [cif, cif_a, cif_b]
    else:
        cifs = [cif]


    #### 
    output = []

    for cf in cifs:

        ## Calculate the ligand positions
        centralIon, ligandPositions, ligandCharge, inv, ligandNames, CartM, CentralIonPos = FindPointGroupSymOps(cf, 
                                                                    mag_ion, Zaxis, 
                                                                    Yaxis, crystalImage,NumIonNeighbors,
                                                                    CoordinationNumber, MaxDistance)
        #print(ligandNames)
        if centralIon in ION_NUMS_RARE_EARTH: # It's a rare earth ion
            if LS_Coupling:
                Lig = LS_Ligands(ion=centralIon, ionPos = [0,0,0], ligandPos = ligandPositions, 
                            SpinOrbitCoupling=LS_Coupling)

            else:
                Lig = Ligands(ion=centralIon, ionPos = [0,0,0], ligandPos = ligandPositions)
            # Create a point charge model, assuming that a mirror plane has been found.
            print('   Creating a point charge model...')
            if ForceImaginary:
                PCM = Lig.PointChargeModel(printB = True, LigandCharge=ligandCharge, suppressminusm = False)
            else:
                PCM = Lig.PointChargeModel(printB = True, LigandCharge=ligandCharge, suppressminusm = inv)

        else: # It's not a rare earth!
            checkTMexist(centralIon)
            if (ionL == None) | (ionS == None):
                raise TypeError('\tplease specify the ionL and ionS values in the importCIF function for '+ centralIon)

            if LS_Coupling: # User-provided SOC
                Lig = LS_Ligands(ion=[centralIon, ionS, ionL], ionPos = [0,0,0], 
                        ligandPos = ligandPositions,  SpinOrbitCoupling=LS_Coupling)
            else: # Look up SOC in a table
                print('    No SOC provided, assuming SOC =', np.around(SPIN_ORBIT_COUPLING_CONSTANTS[centralIon],2), 'meV for '+
                       centralIon +"\n           (if you'd like to adjust this, use the 'LS_Coupling' command).\n")
                
                Lig = LS_Ligands(ion=[centralIon, ionS, ionL], ionPos = [0,0,0], 
                        ligandPos = ligandPositions,  SpinOrbitCoupling=SPIN_ORBIT_COUPLING_CONSTANTS[centralIon])

            if ForceImaginary:
                PCM = Lig.TMPointChargeModel(printB = True, LigandCharge=ligandCharge, suppressminusm = False)
            else:
                PCM = Lig.TMPointChargeModel(printB = True, LigandCharge=ligandCharge, suppressminusm = inv)


        Lig.LatticeTransformM = CartM ## Transformation from lattice ABC to the axes for the bonds
        Lig.LigandNames = ligandNames
        Lig.CentralIonPos = CentralIonPos
        output.append([Lig, PCM])

    if len(output) == 1:
        return Lig, PCM
    else: 
        print('WARNING: more than one ligand position given...\n  '+
            ' outputting [[Ligands1, CFLevels1], [Ligands2, CFLevels2], [Ligands3, CFLevels3]]')
        return output


def checkTMexist(ion):
    """
    Verify transition metal ion has database entries.
    
    Prints warnings if radial integrals, spin-orbit coupling, or
    shell filling data is missing for the specified ion.
    
    Args:
        ion: Ion symbol (e.g., 'Ni2+', 'Fe3+')
    
    Note:
        Non-fatal warnings. User can provide missing data via function arguments.
    """
    if ion not in RADIAL_INTEGRALS_TRANS_METAL:
        print(ion, 'radial integrals are not known by PyCrystalField.')
    if ion not in SPIN_ORBIT_COUPLING_CONSTANTS:
        print(ion, 'spin orbit coupling constant is not known by PyCrystalField, but can be specified with `LS_Coupling`.')
    if ion not in ION_HALF_FILLED and ion not in ION_NOT_HALF_FILLED:
        print(ion, 'shell filling not known by PyCrystalField.')
