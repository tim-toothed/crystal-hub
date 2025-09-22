# Code for computing crystal electric fields

from constants import Jion, JionTM, LambdaConstants # Словари с константами: ионы Rare Earth, ионы Transition Metal, Wybourne-Stevens
from cf_levels import CFLevels, LS_CFLevels
from ligands import Ligands, LS_Ligands 
from import_CIF import importCIF

from thermo_functions import partition_func, Cp_from_CEF # Не используется нигде
from Wybourne_Stevens import WybourneToStevens, StevensToWybourne # Не используется нигде
from rescale_CEF import rescaleCEF # Не используется нигде


print(' '+'*'*55 + '\n'+
     ' *                PyCrystalField 2.3.11               *\n' +
    ' *  Please cite  J. Appl. Cryst. (2021). 54, 356-362   * \n' +
    ' *    <https://doi.org/10.1107/S160057672001554X>      *\n ' + '*'*55+'\n')

