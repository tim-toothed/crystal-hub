import numpy as np
from copy import deepcopy

from .cif_symmetry_import import FindPointGroupSymOps
from .cif_file import CifFile
from .constants import ION_NUMS_RARE_EARTH, SPIN_ORBIT_COUPLING_CONSTANTS, RADIAL_INTEGRALS_TRANS_METAL, ION_HALF_FILLED, ION_NOT_HALF_FILLED
from .ligands import Ligands, LS_Ligands


def importCIF(ciffile, mag_ion=None, Zaxis=None, Yaxis=None, LS_Coupling=None,
              crystalImage=False, NumIonNeighbors=1, ForceImaginary=False,
              ionL=None, ionS=None, CoordinationNumber=None, MaxDistance=None):
    """
    Generate a PyCrystalField point charge model from a CIF file.
    
    Args:
        ciffile: Path to CIF file
        mag_ion: Magnetic ion to analyze (if None, auto-detect first rare earth)
        Zaxis: Optional preferred z-axis
        Yaxis: Optional preferred y-axis
        LS_Coupling: Optional spin-orbit coupling value
        crystalImage: Whether to plot crystal structure
        NumIonNeighbors: Number of nearest neighbor types
        ForceImaginary: Force complex eigenstates
        ionL: Orbital angular momentum (for transition metals)
        ionS: Spin quantum number (for transition metals)
        CoordinationNumber: Optional coordination number
        MaxDistance: Optional maximum distance for neighbors
        
    Returns:
        tuple or list: (Ligands, CFLevels) or list of multiple configurations
    """
    cif = CifFile(ciffile)
    if mag_ion == None:  # take the first rare earth in the cif file as the central ion
        for asuc in cif.asymunitcell:
            if asuc[1].strip('3+') in ['Sm', 'Pm', 'Nd', 'Ce', 'Dy', 'Ho', 'Tm', 'Pr', 'Er', 'Tb', 'Yb']:
                mag_ion = asuc[0]
                print('No mag_ion ion listed, assuming', mag_ion, 'is the central ion.')
                break

    ## Check for multiply defined atoms
    differentPositionsA = []
    differentPositionsB = []
    for ii, at in enumerate(cif.unitcell):
        if at[4] < 0:
            print('negative atom!', ii, at)

        if at[0][-1] in ["'", "B", "b"]:
            differentPositionsA.append(at[0])
            differentPositionsB.append(at[0].replace("'", "").replace("B", "A").replace("b", "a"))

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
        centralIon, ligandPositions, ligandCharge, inv, ligandNames, CartM, CentralIonPos = FindPointGroupSymOps(
            cf, mag_ion, Zaxis, Yaxis, crystalImage, NumIonNeighbors, CoordinationNumber, MaxDistance)
        
        if centralIon in ION_NUMS_RARE_EARTH:  # It's a rare earth ion
            if LS_Coupling:
                Lig = LS_Ligands(ion=centralIon, ionPos=[0, 0, 0], ligandPos=ligandPositions,
                                SpinOrbitCoupling=LS_Coupling)
            else:
                Lig = Ligands(ion=centralIon, ionPos=[0, 0, 0], ligandPos=ligandPositions)
            
            # Create a point charge model, assuming that a mirror plane has been found.
            print('   Creating a point charge model...')
            if ForceImaginary:
                PCM = Lig.PointChargeModel(printB=True, LigandCharge=ligandCharge, suppressminusm=False)
            else:
                PCM = Lig.PointChargeModel(printB=True, LigandCharge=ligandCharge, suppressminusm=inv)

        else:  # It's not a rare earth!
            checkTMexist(centralIon)
            if (ionL == None) | (ionS == None):
                raise TypeError('\tplease specify the ionL and ionS values in the importCIF function for ' + centralIon)

            if LS_Coupling:  # User-provided SOC
                Lig = LS_Ligands(ion=[centralIon, ionS, ionL], ionPos=[0, 0, 0],
                                ligandPos=ligandPositions, SpinOrbitCoupling=LS_Coupling)
            else:  # Look up SOC in a table
                print('    No SOC provided, assuming SOC =', np.around(SPIN_ORBIT_COUPLING_CONSTANTS[centralIon], 2),
                      'meV for ' + centralIon + "\n           (if you'd like to adjust this, use the 'LS_Coupling' command).\n")

                Lig = LS_Ligands(ion=[centralIon, ionS, ionL], ionPos=[0, 0, 0],
                                ligandPos=ligandPositions, SpinOrbitCoupling=SPIN_ORBIT_COUPLING_CONSTANTS[centralIon])

            if ForceImaginary:
                PCM = Lig.TMPointChargeModel(printB=True, LigandCharge=ligandCharge, suppressminusm=False)
            else:
                PCM = Lig.TMPointChargeModel(printB=True, LigandCharge=ligandCharge, suppressminusm=inv)

        Lig.LatticeTransformM = CartM  ## Transformation from lattice ABC to the axes for the bonds
        Lig.LigandNames = ligandNames
        Lig.CentralIonPos = CentralIonPos
        output.append([Lig, PCM])

    if len(output) == 1:
        return Lig, PCM
    else:
        print('WARNING: more than one ligand position given...\n  ' +
              ' outputting [[Ligands1, CFLevels1], [Ligands2, CFLevels2], [Ligands3, CFLevels3]]')
        return output


def checkTMexist(ion):
    """Check if transition metal ion data exists in the database."""
    if ion not in RADIAL_INTEGRALS_TRANS_METAL:
        print(ion, 'radial integrals are not known by PyCrystalField.')
    if ion not in SPIN_ORBIT_COUPLING_CONSTANTS:
        print(ion, 'spin orbit coupling constant is not known by PyCrystalField, but can be specified with `LS_Coupling`.')
    if ion not in ION_HALF_FILLED and ion not in ION_NOT_HALF_FILLED:
        print(ion, 'shell filling not known by PyCrystalField.')
