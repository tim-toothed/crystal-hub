"""
PyCrystalField Package

A package for crystal field calculations in rare earth and transition metal systems.
"""

__version__ = "0.1.0-alpha"
__author__ = "Your Name"

from .cif_file import CifFile
from .ligands import Ligands, LS_Ligands
from .cf_levels import CFLevels, LS_CFLevels
from .import_CIF import importCIF

__all__ = [
    'CifFile',
    'Ligands',
    'LS_Ligands',
    'CFLevels',
    'LS_CFLevels',
    'importCIF',
]