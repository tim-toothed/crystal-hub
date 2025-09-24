# Отредактированные модули (pcf_lib + PyCrystalField)<br>
<br>
1. [wybourne_stevens.py](wybourne_stevens.py)<br>
    **Contents**:<br>
        - `WybourneToStevens` (function)<br>
        - `StevensToWybourne` (function)<br>
    **Inner Dependencies**:<br>
        - `from constants import LambdaConstants, LStheta, theta`<br>
    **Outer Dependencies**:<br>
        - None<br>
<br>
2. [thermo_functions.py](thermo_functions.py)<br>
    **Contents**:<br>
        - `partition_func` (function)<br>
        - `Cp_from_CEF` (function)<br>
    **Inner Dependencies**:<br>
        - `from constants import k_b`<br>
    **Outer Dependencies**:<br>
        - `import numpy as np`<br>
<br>
3. [stevens_operators.py](stevens_operators.py)<br>
    **Contents**:<br>
        - `StevensOp` (function)<br>
        - `LS_StevensOp` (function)<br>
    **Inner Dependencies**:<br>
        - `from operators import Operator`<br>
    **Outer Dependencies**:<br>
        - `import numpy as np`<br>
<br>
4. [rescale_CEF.py](rescale_CEF.py)<br>
    **Contents**:<br>
        - `rescaleCEF` (function)<br>
    **Inner Dependencies**:<br>
        - `from constants import theta, RadialIntegral`<br>
    **Outer Dependencies**:<br>
        - None<br>
<br>
5. [PyCrystalField.py](PyCrystalField.py)<br>
    Main file for all imports<br>
<br>
6. [plot_ligands.py](plot_ligands.py)<br>
    **Contents**:<br>
        - `plotPCF` (function)<br>
        - `exportLigandCif` (function)<br>
        - `atomplot` (class)<br>
            - `__init__` (method)<br>
            - `plotatoms` (method)<br>
            - `plotaxes` (method)<br>
            - `plotabc` (method)<br>
            - `_flatten` (method)<br>
    **Inner Dependencies**:<br>
        - None<br>
    **Outer Dependencies**:<br>
        - `import numpy as np`<br>
        - `import matplotlib.pyplot as plt`<br>
<br>
7. [operators.py](operators.py)<br>
    **Contents**:<br>
        - `Ket` (class)<br>
            - `__init__` (method)<br>
            - `Jz` (method)<br>
            - `Jplus` (method)<br>
            - `Jminus` (method)<br>
            - `Jx` (method)<br>
            - `Jy` (method)<br>
            - `R` (method)<br>
            - `_Rz` (method)<br>
            - `_Ry` (method)<br>
            - `_WignersFormula` (method)<br>
            - `__mul__` (method)<br>
            - `__add__` (method)<br>
        - `Operator` (class)<br>
            - `__init__` (method)<br>
            - `Jz` (method)<br>
            - `Jplus` (method)<br>
            - `Jminus` (method)<br>
            - `Jx` (method)<br>
            - `Jy` (method)<br>
            - `__add__` (method)<br>
            - `__radd__` (method)<br>
            - `__sub__` (method)<br>
            - `__mul__` (method)<br>
            - `__rmul__` (method)<br>
            - `__pow__` (method)<br>
            - `__neg__` (method)<br>
            - `__repr__` (method)<br>
        - `LSOperator` (class)<br>
            - `__init__` (method)<br>
            - `Lz` (method)<br>
            - `Lplus` (method)<br>
            - `Lminus` (method)<br>
            - `Lx` (method)<br>
            - `Ly` (method)<br>
            - `Sz` (method)<br>
            - `Splus` (method)<br>
            - `Sminus` (method)<br>
            - `Sx` (method)<br>
            - `Sy` (method)<br>
            - `__add__` (method)<br>
            - `__radd__` (method)<br>
            - `__sub__` (method)<br>
            - `__mul__` (method)<br>
            - `__rmul__` (method)<br>
            - `__pow__` (method)<br>
            - `__neg__` (method)<br>
            - `__repr__` (method)<br>
            - `magnetization` (method)<br>
            - `susceptibility` (method)<br>
    **Inner Dependencies**:<br>
        - None<br>
    **Outer Dependencies**:<br>
        - `import numpy as np`<br>
<br>
8. [moments_of_inertia.py](moments_of_inertia.py)<br>
    **Contents**:<br>
        - `MomIntertia` (function)<br>
        - `selectZaxisMI` (function)<br>
        - `ContinuousShapeMeasure` (function)<br>
        - `anglesToVector` (function)<br>
        - `rotationMatrix` (function)<br>
        - `rotateArbAxis` (function)<br>
        - `findZaxis_SOM_rotation` (function)<br>
        - `findZaxis` (function)<br>
    **Inner Dependencies**:<br>
        - None<br>
    **Outer Dependencies**:<br>
        - `import numpy as np`<br>
        - `from numba import njit`<br>
        - `from scipy.optimize import minimize`<br>
<br>
9. [ligands.py](ligands.py)<br>
    **Contents**:<br>
        - `Ligands` (class)<br>
            - `__init__` (method)<br>
            - `rotateLigands` (method)<br>
            - `rotateLigandsZ` (method)<br>
            - `_rotateMatrix` (method)<br>
            - `exportCif` (method)<br>
            - `PointChargeModel` (method)<br>
            - `FitChargesNeutrons` (method)<br>
            - `FitCharges` (method)<br>
        - `LS_Ligands` (class)<br>
            - `__init__` (method)<br>
            - `rotateLigands` (method)<br>
            - `rotateLigandsZ` (method)<br>
            - `_rotateMatrix` (method)<br>
            - `PointChargeModel` (method)<br>
            - `TMPointChargeModel` (method)<br>
            - `UnknownTMPointChargeModel` (method)<br>
            - `FitChargesNeutrons` (method)<br>
    **Inner Dependencies**:<br>
        - `from lattice_class import Lattice`<br>
        - `from plot_ligands import exportLigandCif`<br>
        - `from constants import TessHarm, theta, RadialIntegral, Constant, LStheta, PFalpha, PFbeta, RadialIntegral_TM, Jion`<br>
        - `from half_filled import IsHalfFilled`<br>
        - `from stevens_operators import StevensOp, LS_StevensOp`<br>
        - `from create_fit_function import makeFitFunction`<br>
        - `from cf_levels import CFLevels, LS_CFLevels`<br>
        - `from operators import LSOperator`<br>
    **Outer Dependencies**:<br>
        - `import numpy as np`<br>
        - `from scipy import optimize`<br>
<br>
10. [lattice_class.py](lattice_class.py)<br>
    **Contents**:<br>
        - `Lattice` (class)<br>
            - `__init__` (method)<br>
            - `reciplatt` (method)<br>
            - `cartesian` (method)<br>
            - `ABC` (method)<br>
            - `inverseA` (method)<br>
    **Inner Dependencies**:<br>
        - None<br>
    **Outer Dependencies**:<br>
        - `import numpy as np`<br>
<br>
11. [latex_cef_print.py](latex_cef_print.py)<br>
    **Contents**:<br>
        - `printLaTexCEFparams` (function)<br>
    **Inner Dependencies**:<br>
        - None<br>
    **Outer Dependencies**:<br>
        - `import numpy as np`<br>
<br>
12. [import_CIF.py](import_CIF.py)<br>
    **Contents**:<br>
        - `importCIF` (function)<br>
        - `checkTMexist` (function)<br>
    **Inner Dependencies**:<br>
        - `from cifsymmetryimport import FindPointGroupSymOps`<br>
        - `from cif_file import CifFile`<br>
        - `from constants import Jion, SpOrbCoup, TMradialI, HalfList, notHalfList`<br>
        - `from ligands import Ligands, LS_Ligands`<br>
    **Outer Dependencies**:<br>
        - `import numpy as np`<br>
        - `from copy import deepcopy`<br>
<br>
13. [half_filled.py](half_filled.py)<br>
    **Contents**:<br>
        - `IsHalfFilled` (function) <- копирует `notHalfList` (dict) и `HalfList` (dict), которые есть в constants - перенести<br>
    **Inner Dependencies**:<br>
        - None<br>
    **Outer Dependencies**:<br>
        - None<br>
<br>
14. [form_factors.py](form_factors.py)<br>
    **Contents**:<br>
        - `importRE_FF` (function)<br>
        - `RE_FormFactor` (function)<br>
    **Inner Dependencies**:<br>
        - `from constants import Jion`<br>
        - [RE_formfactors.pck] <- перенести в локальное окружение<br>
    **Outer Dependencies**:<br>
        - `import numpy as np`<br>
<br>
15. [create_fit_function.py](create_fit_function.py) <- ERROR: inconsistent use of tabs and spaces in indentation (Беды_с_башкой.py)<br>
    **Contents**:<br>
        - `makeFitFunction` (function)<br>
        - `makeCurveFitFunction` (function)<br>
    **Inner Dependencies**:<br>
        - None<br>
    **Outer Dependencies**:<br>
        - None<br>
<br>
16. [constants.py](constants.py)<br>
    **Contents**:<br>
        - `JionTM` (dict)<br>
        - `Jion` (dict)<br>
        - `LambdaConstants` (dict)      <br>
        - `TESSERAL_CONSTANTS` (dict) <br>
        - `Constant` (function) <- обращается к  TESSERAL_CONSTANTS в global env<br>
        - `TessHarm` (function) <- обращается к  TESSERAL_CONSTANTS в global env  <br>
        - `_tesseral_dispatch` (function) <- обращается к  TESSERAL_CONSTANTS в global env  <br>
        - `PFalpha` (function) <- не используется<br>
        - `PFbeta` (function) <- не используется<br>
        - `PFgamma` (function) <- не используется<br>
        - `LStheta` (function) <- содержит LSThet (dict)<br>
        - `theta` (function) <- содержит Thet (dict)<br>
        - `SPIN_ORBIT_COUPLING_CM` (dict) <- не используется<br>
        - `SpOrbCoup` (dict)<br>
        - `radialI` (dict)<br>
        - `TMradialI` (dict)<br>
        - `RadialIntegral` (function) <- обращается к `radialI` (dict) в global env<br>
        - `RadialIntegral_TM` (function) <- обращается к `TMradialI` (dict) в global env<br>
        - `HalfList` (dict)<br>
        - `notHalfList` (dict)<br>
        - `ahc` (float)<br>
        - `a0` (float)<br>
        - `muB` (float)<br>
        - `k_B` (float)<br>
    **Inner Dependencies**:<br>
        - ``<br>
    **Outer Dependencies**:<br>
        - `import numpy as np`<br>
        - `from typing import Union, Dict`<br>
<br>
17. [cif_symmetry_import.py](cif_symmetry_import.py) <- ERROR: inconsistent use of tabs and spaces in indentation<br>
    **Contents**:<br>
        - `FindPointGroupSymOps` (function)<br>
        - `findRotationAxis` (function)<br>
        - `makeSymOpMatrix` (function)<br>
    **Inner Dependencies**:<br>
        - `from plot_ligands import plotPCF`<br>
        - `from moments_of_inertia import findZaxis`<br>
    **Outer Dependencies**:<br>
        - `import numpy as np`<br>
<br>
18. [cif_file.py](cif_file.py)<br>
    **Contents**:<br>
        - `CifFile` (Class)<br>
            - `__init__` (method)<br>
            - `SymOperate` (method)<br>
            - `MakeUnitCell` (method)<br>
            - `StructureFactor` (method)<br>
            - `MultipleScattering` (method)<br>
            - `_destringify` (method)<br>
            - `_defractionify` (method)<br>
            - `_duplicaterow` (method)<br>
            - `_NumElements` (method)<br>
            - `_kvector` (method)<br>
    **Inner Dependencies**:<br>
        - `from lattice_class import lattice`<br>
    **Outer Dependencies**:<br>
        - `import numpy as np`<br>
        - `from copy import deepcopy`<br>
<br>
19. [cf_levels.py](cf_levels.py)<br>
    **Contents**:<br>
        - `CFLevels` (class)<br>
            - `__init__` (method)<br>
            - `Bdict` (method)<br>
            - `Hamiltonian` (method)<br>
            - `newCoeff` (method)<br>
            - `diagonalize` (method)<br>
            - `diagonalize_banded` (method)<br>
            - `_findbands` (method)<br>
            - `transitionIntensity` (method)<br>
            - `neutronSpectrum` (method)<br>
            - `neutronSpectrum_customLineshape` (method)<br>
            - `normalizedNeutronSpectrum` (method)<br>
            - `normalizedNeutronSpectrum_customLineshape` (method)<br>
            - `neutronSpectrum2D` (method)<br>
            - `normalizedNeutronSpectrum2D` (method)<br>
            - `_transition` (method)<br>
            - `_lorentzian` (method)<br>
            - `_voigt` (method)<br>
            - `_Re` (method)<br>
            - `printEigenvectors` (method)<br>
            - `printLaTexEigenvectors` (method)<br>
            - `gsExpectation` (method)<br>
            - `magnetization` (method)<br>
            - `susceptibility` (method)<br>
            - `susceptibilityPert` (method)<br>
            - `gtensor` (method)<br>
            - `gtensorzeeman` (method)<br>
            - `fitdata` (method)<br>
            - `fitdata_GlobalOpt` (method)<br>
            - `testEigenvectors` (method)<br>
        - `_spec` (list) <- ???<br>
        - `OpticalTransition` (class)<br>
            - `__init__` (method)<br>
            - `transition_strength` (method)               <br>
        - `LandeGFactor` (function)<br>
        - `LS_CFLevels` (class)<br>
            - `__init__` (method)<br>
            - `Bdict` (method)<br>
            - `Hamiltonian` (method)<br>
            - `newCoeff` (method)<br>
            - `diagonalize` (method)<br>
            - `_findbands` (method)<br>
            - `neutronSpectrum` (method)<br>
            - `neutronSpectrum2D` (method)<br>
            - `normalizedNeutronSpectrum` (method)<br>
            - `_transition` (method)<br>
            - `_lorentzian` (method)<br>
            - `_voigt` (method)<br>
            - `_Re` (method)<br>
            - `printEigenvectors` (method)<br>
            - `gsExpectation` (method)<br>
            - `magnetization` (method)<br>
            - `susceptibility` (method)<br>
            - `susceptibilityDeriv` (method)<br>
            - `magnetizationDeriv` (method)<br>
            - `gtensor` (method)<br>
            - `gtensorperturb` (method)<br>
            - `fitdata` (method)<br>
            - `printLaTexEigenvectors` (method)<br>
    **Inner Dependencies**:<br>
        - `from constants import Jion`<br>
        - `from form_factors import RE_FormFactor`<br>
        - `from create_fit_function import makeFitFunction`<br>
        - `from operators import Ket, Operator, LSOperator`<br>
        - `from stevens_operators import StevensOp, LS_StevensOp`<br>
    **Outer Dependencies**:<br>
        - `import numpy as np`<br>
        - `from numba import jitclass, float64`<br>
        - `from scipy import optimize`<br>
        - `import scipy.linalg as LA`<br>
        - `from scipy.special import wofz`