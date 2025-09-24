# Отредактированные модули (pcf_lib + PyCrystalField)

<details>
<summary>1. wybourne_stevens.py</summary>

**Contents**:
- `WybourneToStevens` (function)
- `StevensToWybourne` (function)

**Inner Dependencies**:
- `from constants import LambdaConstants, LStheta, theta`

**Outer Dependencies**:
- None
</details>

<details>
<summary>2. thermo_functions.py</summary>

**Contents**:
- `partition_func` (function)
- `Cp_from_CEF` (function)

**Inner Dependencies**:
- `from constants import k_b`

**Outer Dependencies**:
- `import numpy as np`
</details>

<details>
<summary>3. stevens_operators.py</summary>

**Contents**:
- `StevensOp` (function)
- `LS_StevensOp` (function)

**Inner Dependencies**:
- `from operators import Operator`

**Outer Dependencies**:
- `import numpy as np`
</details>

<details>
<summary>4. rescale_CEF.py</summary>

**Contents**:
- `rescaleCEF` (function)

**Inner Dependencies**:
- `from constants import theta, RadialIntegral`

**Outer Dependencies**:
- None
</details>

<details>
<summary>5. PyCrystalField.py</summary>

Main file for all imports
</details>

<details>
<summary>6. plot_ligands.py</summary>

**Contents**:
- `plotPCF` (function)
- `exportLigandCif` (function)
- `atomplot` (class)
  - `__init__` (method)
  - `plotatoms` (method)
  - `plotaxes` (method)
  - `plotabc` (method)
  - `_flatten` (method)

**Inner Dependencies**:
- None

**Outer Dependencies**:
- `import numpy as np`
- `import matplotlib.pyplot as plt`
</details>

<details>
<summary>7. operators.py</summary>

**Contents**:
- `Ket` (class)
  - `__init__` (method)
  - `Jz` (method)
  - `Jplus` (method)
  - `Jminus` (method)
  - `Jx` (method)
  - `Jy` (method)
  - `R` (method)
  - `_Rz` (method)
  - `_Ry` (method)
  - `_WignersFormula` (method)
  - `__mul__` (method)
  - `__add__` (method)
- `Operator` (class)
  - `__init__` (method)
  - `Jz` (method)
  - `Jplus` (method)
  - `Jminus` (method)
  - `Jx` (method)
  - `Jy` (method)
  - `__add__` (method)
  - `__radd__` (method)
  - `__sub__` (method)
  - `__mul__` (method)
  - `__rmul__` (method)
  - `__pow__` (method)
  - `__neg__` (method)
  - `__repr__` (method)
- `LSOperator` (class)
  - `__init__` (method)
  - `Lz` (method)
  - `Lplus` (method)
  - `Lminus` (method)
  - `Lx` (method)
  - `Ly` (method)
  - `Sz` (method)
  - `Splus` (method)
  - `Sminus` (method)
  - `Sx` (method)
  - `Sy` (method)
  - `__add__` (method)
  - `__radd__` (method)
  - `__sub__` (method)
  - `__mul__` (method)
  - `__rmul__` (method)
  - `__pow__` (method)
  - `__neg__` (method)
  - `__repr__` (method)
  - `magnetization` (method)
  - `susceptibility` (method)

**Inner Dependencies**:
- None

**Outer Dependencies**:
- `import numpy as np`
</details>

<details>
<summary>8. moments_of_inertia.py</summary>

**Contents**:
- `MomIntertia` (function)
- `selectZaxisMI` (function)
- `ContinuousShapeMeasure` (function)
- `anglesToVector` (function)
- `rotationMatrix` (function)
- `rotateArbAxis` (function)
- `findZaxis_SOM_rotation` (function)
- `findZaxis` (function)

**Inner Dependencies**:
- None

**Outer Dependencies**:
- `import numpy as np`
- `from numba import njit`
- `from scipy.optimize import minimize`
</details>

<details>
<summary>9. ligands.py</summary>

**Contents**:
- `Ligands` (class)
  - `__init__` (method)
  - `rotateLigands` (method)
  - `rotateLigandsZ` (method)
  - `_rotateMatrix` (method)
  - `exportCif` (method)
  - `PointChargeModel` (method)
  - `FitChargesNeutrons` (method)
  - `FitCharges` (method)
- `LS_Ligands` (class)
  - `__init__` (method)
  - `rotateLigands` (method)
  - `rotateLigandsZ` (method)
  - `_rotateMatrix` (method)
  - `PointChargeModel` (method)
  - `TMPointChargeModel` (method)
  - `UnknownTMPointChargeModel` (method)
  - `FitChargesNeutrons` (method)

**Inner Dependencies**:
- `from lattice_class import Lattice`
- `from plot_ligands import exportLigandCif`
- `from constants import TessHarm, theta, RadialIntegral, Constant, LStheta, PFalpha, PFbeta, RadialIntegral_TM, Jion`
- `from half_filled import IsHalfFilled`
- `from stevens_operators import StevensOp, LS_StevensOp`
- `from create_fit_function import makeFitFunction`
- `from cf_levels import CFLevels, LS_CFLevels`
- `from operators import LSOperator`

**Outer Dependencies**:
- `import numpy as np`
- `from scipy import optimize`
</details>

<details>
<summary>10. lattice_class.py</summary>

**Contents**:
- `Lattice` (class)
  - `__init__` (method)
  - `reciplatt` (method)
  - `cartesian` (method)
  - `ABC` (method)
  - `inverseA` (method)

**Inner Dependencies**:
- None

**Outer Dependencies**:
- `import numpy as np`
</details>

<details>
<summary>11. latex_cef_print.py</summary>

**Contents**:
- `printLaTexCEFparams` (function)

**Inner Dependencies**:
- None

**Outer Dependencies**:
- `import numpy as np`
</details>

<details>
<summary>12. import_CIF.py</summary>

**Contents**:
- `importCIF` (function)
- `checkTMexist` (function)

**Inner Dependencies**:
- `from cifsymmetryimport import FindPointGroupSymOps`
- `from cif_file import CifFile`
- `from constants import Jion, SpOrbCoup, TMradialI, HalfList, notHalfList`
- `from ligands import Ligands, LS_Ligands`

**Outer Dependencies**:
- `import numpy as np`
- `from copy import deepcopy`
</details>

<details>
<summary>13. half_filled.py</summary>

**Contents**:
- `IsHalfFilled` (function) <- копирует `notHalfList` (dict) и `HalfList` (dict), которые есть в constants - перенести

**Inner Dependencies**:
- None

**Outer Dependencies**:
- None
</details>

<details>
<summary>14. form_factors.py</summary>

**Contents**:
- `importRE_FF` (function)
- `RE_FormFactor` (function)

**Inner Dependencies**:
- `from constants import Jion`
- [RE_formfactors.pck] <- перенести в локальное окружение

**Outer Dependencies**:
- `import numpy as np`
</details>

<details>
<summary>15. create_fit_function.py</summary>

**ERROR**: inconsistent use of tabs and spaces in indentation (Беды_с_башкой.py)

**Contents**:
- `makeFitFunction` (function)
- `makeCurveFitFunction` (function)

**Inner Dependencies**:
- None

**Outer Dependencies**:
- None
</details>

<details>
<summary>16. constants.py</summary>

**Contents**:
- `JionTM` (dict)
- `Jion` (dict)
- `LambdaConstants` (dict)
- `TESSERAL_CONSTANTS` (dict)
- `Constant` (function) <- обращается к  TESSERAL_CONSTANTS в global env
- `TessHarm` (function) <- обращается к  TESSERAL_CONSTANTS в global env
- `_tesseral_dispatch` (function) <- обращается к  TESSERAL_CONSTANTS в global env
- `PFalpha` (function) <- не используется
- `PFbeta` (function) <- не используется
- `PFgamma` (function) <- не используется
- `LStheta` (function) <- содержит LSThet (dict)
- `theta` (function) <- содержит Thet (dict)
- `SPIN_ORBIT_COUPLING_CM` (dict) <- не используется
- `SpOrbCoup` (dict)
- `radialI` (dict)
- `TMradialI` (dict)
- `RadialIntegral` (function) <- обращается к `radialI` (dict) в global env
- `RadialIntegral_TM` (function) <- обращается к `TMradialI` (dict) в global env
- `HalfList` (dict)
- `notHalfList` (dict)
- `ahc` (float)
- `a0` (float)
- `muB` (float)
- `k_B` (float)

**Inner Dependencies**:
- ``

**Outer Dependencies**:
- `import numpy as np`
- `from typing import Union, Dict`
</details>

<details>
<summary>17. cif_symmetry_import.py</summary>

**ERROR**: inconsistent use of tabs and spaces in indentation

**Contents**:
- `FindPointGroupSymOps` (function)
- `findRotationAxis` (function)
- `makeSymOpMatrix` (function)

**Inner Dependencies**:
- `from plot_ligands import plotPCF`
- `from moments_of_inertia import findZaxis`

**Outer Dependencies**:
- `import numpy as np`
</details>

<details>
<summary>18. cif_file.py</summary>

**Contents**:
- `CifFile` (Class)
  - `__init__` (method)
  - `SymOperate` (method)
  - `MakeUnitCell` (method)
  - `StructureFactor` (method)
  - `MultipleScattering` (method)
  - `_destringify` (method)
  - `_defractionify` (method)
  - `_duplicaterow` (method)
  - `_NumElements` (method)
  - `_kvector` (method)

**Inner Dependencies**:
- `from lattice_class import lattice`

**Outer Dependencies**:
- `import numpy as np`
- `from copy import deepcopy`
</details>

<details>
<summary>19. cf_levels.py</summary>

**Contents**:
- `CFLevels` (class)
  - `__init__` (method)
  - `Bdict` (method)
  - `Hamiltonian` (method)
  - `newCoeff` (method)
  - `diagonalize` (method)
  - `diagonalize_banded` (method)
  - `_findbands` (method)
  - `transitionIntensity` (method)
  - `neutronSpectrum` (method)
  - `neutronSpectrum_customLineshape` (method)
  - `normalizedNeutronSpectrum` (method)
  - `normalizedNeutronSpectrum_customLineshape` (method)
  - `neutronSpectrum2D` (method)
  - `normalizedNeutronSpectrum2D` (method)
  - `_transition` (method)
  - `_lorentzian` (method)
  - `_voigt` (method)
  - `_Re` (method)
  - `printEigenvectors` (method)
  - `printLaTexEigenvectors` (method)
  - `gsExpectation` (method)
  - `magnetization` (method)
  - `susceptibility` (method)
  - `susceptibilityPert` (method)
  - `gtensor` (method)
  - `gtensorzeeman` (method)
  - `fitdata` (method)
  - `fitdata_GlobalOpt` (method)
  - `testEigenvectors` (method)
- `_spec` (list) <- ???
- `OpticalTransition` (class)
  - `__init__` (method)
  - `transition_strength` (method)
- `LandeGFactor` (function)
- `LS_CFLevels` (class)
  - `__init__` (method)
  - `Bdict` (method)
  - `Hamiltonian` (method)
  - `newCoeff` (method)
  - `diagonalize` (method)
  - `_findbands` (method)
  - `neutronSpectrum` (method)
  - `neutronSpectrum2D` (method)
  - `normalizedNeutronSpectrum` (method)
  - `_transition` (method)
  - `_lorentzian` (method)
  - `_voigt` (method)
  - `_Re` (method)
  - `printEigenvectors` (method)
  - `gsExpectation` (method)
  - `magnetization` (method)
  - `susceptibility` (method)
  - `susceptibilityDeriv` (method)
  - `magnetizationDeriv` (method)
  - `gtensor` (method)
  - `gtensorperturb` (method)
  - `fitdata` (method)
  - `printLaTexEigenvectors` (method)

**Inner Dependencies**:
- `from constants import Jion`
- `from form_factors import RE_FormFactor`
- `from create_fit_function import makeFitFunction`
- `from operators import Ket, Operator, LSOperator`
- `from stevens_operators import StevensOp, LS_StevensOp`

**Outer Dependencies**:
- `import numpy as np`
- `from numba import jitclass, float64`
- `from scipy import optimize`
- `import scipy.linalg as LA`
- `from scipy.special import wofz`
</details>