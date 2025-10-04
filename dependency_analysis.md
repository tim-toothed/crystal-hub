# 📊 Python Package Dependency Analysis Report

**Generated:** 2025-10-05 02:16:21

**Package:** `edited_package`

---

## 📑 Table of Contents

- [cf_levels.py](#cf_levels-py)
- [cif_file.py](#cif_file-py)
- [cif_symmetry_import.py](#cif_symmetry_import-py)
- [constants.py](#constants-py)
- [create_fit_function.py](#create_fit_function-py)
- [form_factors.py](#form_factors-py)
- [fundamental_operators.py](#fundamental_operators-py)
- [import_CIF.py](#import_CIF-py)
- [lattice_class.py](#lattice_class-py)
- [ligands.py](#ligands-py)
- [moments_of_inertia.py](#moments_of_inertia-py)
- [plotPCF.py](#plotPCF-py)
- [rescaleCEF.py](#rescaleCEF-py)
- [stevens_operators.py](#stevens_operators-py)
- [thermo_functions.py](#thermo_functions-py)
- [wybourne_stevens.py](#wybourne_stevens-py)
- [__init__.py](#__init__-py)

---

## 📁 cf_levels.py {#cf_levels-py}

### 🏛️ Classes

#### **CFLevels** *(line 15)*

- **Used by:** 3 files


**Methods:**

- `__init__()` *(line 17)* - 0 usages
- `Bdict()` *(line 29)* - 5 usages
  <details><summary>View all 5 usages</summary>

  - `wybourne_stevens.py`: function:WybourneToStevens (line 31)
  - `wybourne_stevens.py`: function:StevensToWybourne (line 57)
  - `wybourne_stevens.py`: function:_convert_parameters (line 77)
  - `wybourne_stevens.py`: function:_convert_parameters (line 84)
  - `wybourne_stevens.py`: function:_convert_parameters (line 86)
  </details>

- `Hamiltonian()` *(line 45)* - 4 usages
  <details><summary>View all 4 usages</summary>

  - `ligands.py`: class:Ligands → method:PointChargeModel (line 270)
  - `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 607)
  - `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 703)
  - `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 799)
  </details>

- `newCoeff()` *(line 53)* - 0 usages
- `diagonalize()` *(line 58)* - 14 usages
  <details><summary>View all 14 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:newCoeff (line 56)
  - `cf_levels.py`: class:CFLevels → method:neutronSpectrum (line 134)
  - `cf_levels.py`: class:CFLevels → method:neutronSpectrum_customLineshape (line 179)
  - `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum (line 211)
  - `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum_customLineshape (line 242)
  - `cf_levels.py`: class:CFLevels → method:printEigenvectors (line 331)
  - `cf_levels.py`: class:CFLevels → method:printLaTexEigenvectors (line 347)
  - `cf_levels.py`: class:LS_CFLevels → method:newCoeff (line 852)
  - `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum (line 894)
  - `cf_levels.py`: class:LS_CFLevels → method:normalizedNeutronSpectrum (line 943)
  - `cf_levels.py`: class:LS_CFLevels → method:printEigenvectors (line 1027)
  - `cf_levels.py`: class:LS_CFLevels → method:printLaTexEigenvectors (line 1393)
  - `ligands.py`: class:Ligands → method:FitCharges (line 347)
  - `ligands.py`: class:LS_Ligands → method:FitChargesNeutrons (line 844)
  </details>

- `diagonalize_banded()` *(line 80)* - 1 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:gtensor (line 543)
  </details>

- `_findbands()` *(line 98)* - 3 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:diagonalize (line 68)
  - `cf_levels.py`: class:CFLevels → method:diagonalize_banded (line 87)
  - `cf_levels.py`: class:LS_CFLevels → method:diagonalize (line 862)
  </details>

- `transitionIntensity()` *(line 114)* - 0 usages
- `neutronSpectrum()` *(line 126)* - 2 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:neutronSpectrum2D (line 266)
  - `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum2D (line 928)
  </details>

- `neutronSpectrum_customLineshape()` *(line 169)* - 0 usages
- `normalizedNeutronSpectrum()` *(line 203)* - 1 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum2D (line 275)
  </details>

- `normalizedNeutronSpectrum_customLineshape()` *(line 234)* - 0 usages
- `neutronSpectrum2D()` *(line 265)* - 0 usages
- `normalizedNeutronSpectrum2D()` *(line 274)* - 0 usages
- `_transition()` *(line 284)* - 3 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:testEigenvectors (line 699)
  - `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum (line 914)
  - `cf_levels.py`: class:LS_CFLevels → method:normalizedNeutronSpectrum (line 962)
  </details>

- `_lorentzian()` *(line 304)* - 0 usages
- `_voigt()` *(line 307)* - 4 usages
  <details><summary>View all 4 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:neutronSpectrum (line 152)
  - `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum (line 229)
  - `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum (line 917)
  - `cf_levels.py`: class:LS_CFLevels → method:normalizedNeutronSpectrum (line 965)
  </details>

- `_Re()` *(line 314)* - 14 usages
  <details><summary>View all 14 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:printEigenvectors (line 339)
  - `cf_levels.py`: class:CFLevels → method:printEigenvectors (line 339)
  - `cf_levels.py`: class:CFLevels → method:printLaTexEigenvectors (line 367)
  - `cf_levels.py`: class:CFLevels → method:printLaTexEigenvectors (line 368)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 379)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 380)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 381)
  - `cf_levels.py`: class:LS_CFLevels → method:printEigenvectors (line 1035)
  - `cf_levels.py`: class:LS_CFLevels → method:printEigenvectors (line 1035)
  - `cf_levels.py`: class:LS_CFLevels → method:gsExpectation (line 1044)
  - `cf_levels.py`: class:LS_CFLevels → method:gsExpectation (line 1045)
  - `cf_levels.py`: class:LS_CFLevels → method:gsExpectation (line 1046)
  - `cf_levels.py`: class:LS_CFLevels → method:printLaTexEigenvectors (line 1422)
  - `cf_levels.py`: class:LS_CFLevels → method:printLaTexEigenvectors (line 1423)
  </details>

- `printEigenvectors()` *(line 326)* - 0 usages
- `printLaTexEigenvectors()` *(line 342)* - 0 usages
- `gsExpectation()` *(line 372)* - 0 usages
- `magnetization()` *(line 386)* - 48 usages
  <details><summary>View all 48 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 463)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 464)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 465)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 466)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 472)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 473)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 474)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 475)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 481)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 482)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 483)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 484)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 492)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 493)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 494)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 495)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1107)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1108)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1109)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1110)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1116)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1117)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1118)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1119)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1125)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1126)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1127)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1128)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1136)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1137)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1138)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1139)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 653)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 654)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 655)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 656)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 662)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 663)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 664)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 665)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 671)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 672)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 673)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 674)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 682)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 683)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 684)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 685)
  </details>

- `susceptibility()` *(line 449)* - 0 usages
- `susceptibilityPert()` *(line 502)* - 0 usages
- `gtensor()` *(line 540)* - 0 usages
- `gtensorzeeman()` *(line 596)* - 0 usages
- `fitdata()` *(line 643)* - 0 usages
- `fitdata_GlobalOpt()` *(line 667)* - 0 usages
- `testEigenvectors()` *(line 688)* - 0 usages

<details><summary>📍 View all 3 class usages</summary>

- `ligands.py`: import-statement (line 27)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 270)
- `__init__.py`: import-statement (line 12)
</details>

#### **opttransition** *(line 724)*

- **Used by:** 0 file


**Methods:**

- `__init__()` *(line 725)* - 0 usages
- `transition()` *(line 733)* - 0 usages
#### **LS_CFLevels** *(line 768)*

- **Used by:** 5 files


**Methods:**

- `__init__()` *(line 770)* - 0 usages
- `Bdict()` *(line 809)* - 5 usages
  <details><summary>View all 5 usages</summary>

  - `wybourne_stevens.py`: function:WybourneToStevens (line 31)
  - `wybourne_stevens.py`: function:StevensToWybourne (line 57)
  - `wybourne_stevens.py`: function:_convert_parameters (line 77)
  - `wybourne_stevens.py`: function:_convert_parameters (line 84)
  - `wybourne_stevens.py`: function:_convert_parameters (line 86)
  </details>

- `Hamiltonian()` *(line 824)* - 4 usages
  <details><summary>View all 4 usages</summary>

  - `ligands.py`: class:Ligands → method:PointChargeModel (line 270)
  - `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 607)
  - `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 703)
  - `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 799)
  </details>

- `newCoeff()` *(line 849)* - 0 usages
- `diagonalize()` *(line 855)* - 14 usages
  <details><summary>View all 14 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:newCoeff (line 56)
  - `cf_levels.py`: class:CFLevels → method:neutronSpectrum (line 134)
  - `cf_levels.py`: class:CFLevels → method:neutronSpectrum_customLineshape (line 179)
  - `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum (line 211)
  - `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum_customLineshape (line 242)
  - `cf_levels.py`: class:CFLevels → method:printEigenvectors (line 331)
  - `cf_levels.py`: class:CFLevels → method:printLaTexEigenvectors (line 347)
  - `cf_levels.py`: class:LS_CFLevels → method:newCoeff (line 852)
  - `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum (line 894)
  - `cf_levels.py`: class:LS_CFLevels → method:normalizedNeutronSpectrum (line 943)
  - `cf_levels.py`: class:LS_CFLevels → method:printEigenvectors (line 1027)
  - `cf_levels.py`: class:LS_CFLevels → method:printLaTexEigenvectors (line 1393)
  - `ligands.py`: class:Ligands → method:FitCharges (line 347)
  - `ligands.py`: class:LS_Ligands → method:FitChargesNeutrons (line 844)
  </details>

- `_findbands()` *(line 874)* - 3 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:diagonalize (line 68)
  - `cf_levels.py`: class:CFLevels → method:diagonalize_banded (line 87)
  - `cf_levels.py`: class:LS_CFLevels → method:diagonalize (line 862)
  </details>

- `neutronSpectrum()` *(line 889)* - 2 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:neutronSpectrum2D (line 266)
  - `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum2D (line 928)
  </details>

- `neutronSpectrum2D()` *(line 927)* - 0 usages
- `normalizedNeutronSpectrum()` *(line 937)* - 1 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum2D (line 275)
  </details>

- `_transition()` *(line 971)* - 3 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:testEigenvectors (line 699)
  - `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum (line 914)
  - `cf_levels.py`: class:LS_CFLevels → method:normalizedNeutronSpectrum (line 962)
  </details>

- `_lorentzian()` *(line 998)* - 0 usages
- `_voigt()` *(line 1002)* - 4 usages
  <details><summary>View all 4 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:neutronSpectrum (line 152)
  - `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum (line 229)
  - `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum (line 917)
  - `cf_levels.py`: class:LS_CFLevels → method:normalizedNeutronSpectrum (line 965)
  </details>

- `_Re()` *(line 1009)* - 14 usages
  <details><summary>View all 14 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:printEigenvectors (line 339)
  - `cf_levels.py`: class:CFLevels → method:printEigenvectors (line 339)
  - `cf_levels.py`: class:CFLevels → method:printLaTexEigenvectors (line 367)
  - `cf_levels.py`: class:CFLevels → method:printLaTexEigenvectors (line 368)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 379)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 380)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 381)
  - `cf_levels.py`: class:LS_CFLevels → method:printEigenvectors (line 1035)
  - `cf_levels.py`: class:LS_CFLevels → method:printEigenvectors (line 1035)
  - `cf_levels.py`: class:LS_CFLevels → method:gsExpectation (line 1044)
  - `cf_levels.py`: class:LS_CFLevels → method:gsExpectation (line 1045)
  - `cf_levels.py`: class:LS_CFLevels → method:gsExpectation (line 1046)
  - `cf_levels.py`: class:LS_CFLevels → method:printLaTexEigenvectors (line 1422)
  - `cf_levels.py`: class:LS_CFLevels → method:printLaTexEigenvectors (line 1423)
  </details>

- `printEigenvectors()` *(line 1022)* - 0 usages
- `gsExpectation()` *(line 1038)* - 0 usages
- `magnetization()` *(line 1051)* - 48 usages
  <details><summary>View all 48 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 463)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 464)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 465)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 466)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 472)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 473)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 474)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 475)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 481)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 482)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 483)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 484)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 492)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 493)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 494)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 495)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1107)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1108)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1109)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1110)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1116)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1117)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1118)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1119)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1125)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1126)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1127)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1128)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1136)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1137)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1138)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1139)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 653)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 654)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 655)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 656)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 662)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 663)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 664)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 665)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 671)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 672)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 673)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 674)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 682)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 683)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 684)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 685)
  </details>

- `susceptibility()` *(line 1097)* - 0 usages
- `susceptibilityDeriv()` *(line 1146)* - 0 usages
- `magnetizationDeriv()` *(line 1185)* - 12 usages
  <details><summary>View all 12 usages</summary>

  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1156)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1157)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1158)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1159)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1165)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1166)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1167)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1168)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1174)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1175)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1176)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1177)
  </details>

- `gtensor()` *(line 1264)* - 0 usages
- `gtensorperturb()` *(line 1328)* - 0 usages
- `fitdata()` *(line 1359)* - 0 usages
- `printLaTexEigenvectors()` *(line 1388)* - 0 usages

<details><summary>📍 View all 5 class usages</summary>

- `ligands.py`: import-statement (line 27)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 607)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 703)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 799)
- `__init__.py`: import-statement (line 12)
</details>

### 📊 Functions

#### **LandeGFactor()** *(line 743)*

- **Used by:** 0 file

---

## 📁 cif_file.py {#cif_file-py}

### 🏛️ Classes

#### **CifFile** *(line 7)*

- **Used by:** 4 files


**Methods:**

- `__init__()` *(line 9)* - 0 usages
- `SymOperate()` *(line 156)* - 2 usages
  <details><summary>View usages</summary>

  - `cif_file.py`: class:CifFile → method:MakeUnitCell (line 183)
  - `cif_symmetry_import.py`: function:FindPointGroupSymOps (line 48)
  </details>

- `MakeUnitCell()` *(line 178)* - 1 usages
  <details><summary>View usages</summary>

  - `cif_file.py`: class:CifFile → method:__init__ (line 148)
  </details>

- `StructureFactor()` *(line 216)* - 0 usages
- `MultipleScattering()` *(line 258)* - 0 usages
- `_destringify()` *(line 303)* - 11 usages
  <details><summary>View all 11 usages</summary>

  - `cif_file.py`: class:CifFile → method:__init__ (line 30)
  - `cif_file.py`: class:CifFile → method:__init__ (line 32)
  - `cif_file.py`: class:CifFile → method:__init__ (line 34)
  - `cif_file.py`: class:CifFile → method:__init__ (line 36)
  - `cif_file.py`: class:CifFile → method:__init__ (line 38)
  - `cif_file.py`: class:CifFile → method:__init__ (line 40)
  - `cif_file.py`: class:CifFile → method:__init__ (line 84)
  - `cif_file.py`: class:CifFile → method:__init__ (line 85)
  - `cif_file.py`: class:CifFile → method:__init__ (line 86)
  - `cif_file.py`: class:CifFile → method:__init__ (line 88)
  - `cif_file.py`: class:CifFile → method:__init__ (line 90)
  </details>

- `_defractionify()` *(line 321)* - 3 usages
  <details><summary>View usages</summary>

  - `cif_file.py`: class:CifFile → method:SymOperate (line 165)
  - `cif_file.py`: class:CifFile → method:SymOperate (line 166)
  - `cif_file.py`: class:CifFile → method:SymOperate (line 167)
  </details>

- `_duplicaterow()` *(line 325)* - 1 usages
  <details><summary>View usages</summary>

  - `cif_file.py`: class:CifFile → method:MakeUnitCell (line 191)
  </details>

- `_NumElements()` *(line 338)* - 1 usages
  <details><summary>View usages</summary>

  - `cif_file.py`: class:CifFile → method:StructureFactor (line 222)
  </details>

- `_kvector()` *(line 345)* - 1 usages
  <details><summary>View usages</summary>

  - `cif_file.py`: class:CifFile → method:MultipleScattering (line 261)
  </details>


<details><summary>📍 View all 4 class usages</summary>

- `import_CIF.py`: import-statement (line 5)
- `import_CIF.py`: function:importCIF (line 33)
- `import_CIF.py`: function:importCIF (line 33)
- `__init__.py`: import-statement (line 10)
</details>

---

## 📁 cif_symmetry_import.py {#cif_symmetry_import-py}

### 📊 Functions

#### **FindPointGroupSymOps()** *(line 8)*

- **Used by:** 3 files

<details><summary>View all 3 usages</summary>

- `import_CIF.py`: import-statement (line 4)
- `import_CIF.py`: function:importCIF (line 80)
- `import_CIF.py`: function:importCIF (line 80)
</details>

#### **findRotationAxis()** *(line 332)*

- **Used by:** 0 file

#### **makeSymOpMatrix()** *(line 382)*

- **Used by:** 0 file

---

## 📁 constants.py {#constants-py}

### 📊 Functions

#### **Constant()** *(line 100)*

- **Used by:** 9 files

<details><summary>View all 9 usages</summary>

- `ligands.py`: import-statement (line 24)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 258)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 258)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 592)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 592)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 687)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 687)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 784)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 784)
</details>

#### **calculate_tesseral_harmonic()** *(line 105)*

- **Used by:** 9 files

<details><summary>View all 9 usages</summary>

- `ligands.py`: import-statement (line 24)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 254)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 254)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 588)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 588)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 683)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 683)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 780)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 780)
</details>

#### **_tesseral_dispatch()** *(line 117)*

- **Used by:** 0 file

#### **PFalpha()** *(line 187)*

- **Used by:** 5 files

<details><summary>View all 5 usages</summary>

- `ligands.py`: import-statement (line 24)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 671)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 671)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 768)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 768)
</details>

#### **PFbeta()** *(line 207)*

- **Used by:** 5 files

<details><summary>View all 5 usages</summary>

- `ligands.py`: import-statement (line 24)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 672)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 672)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 769)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 769)
</details>

#### **PFgamma()** *(line 238)*

- **Used by:** 0 file

#### **LStheta()** *(line 288)*

- **Used by:** 5 files

<details><summary>View all 5 usages</summary>

- `ligands.py`: import-statement (line 24)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 592)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 592)
- `wybourne_stevens.py`: import-statement (line 1)
- `wybourne_stevens.py`: function:_convert_parameters (line 75)
</details>

#### **theta()** *(line 310)*

- **Used by:** 23 files

<details><summary>View all 23 usages</summary>

- `fundamental_operators.py`: class:Ket → method:_Rz (line 131)
- `ligands.py`: import-statement (line 24)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 258)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 258)
- `moments_of_inertia.py`: function:anglesToVector (line 53)
- `moments_of_inertia.py`: function:anglesToVector (line 53)
- `moments_of_inertia.py`: function:anglesToVector (line 53)
- `moments_of_inertia.py`: function:rotationMatrix (line 58)
- `moments_of_inertia.py`: function:rotateArbAxis (line 73)
- `moments_of_inertia.py`: function:fitfun (line 84)
- `moments_of_inertia.py`: function:fitfun (line 85)
- `plotPCF.py`: class:atomplot → method:__init__ (line 74)
- `plotPCF.py`: class:atomplot → method:__init__ (line 74)
- `plotPCF.py`: class:atomplot → method:__init__ (line 77)
- `plotPCF.py`: class:atomplot → method:__init__ (line 77)
- `plotPCF.py`: class:atomplot → method:__init__ (line 77)
- `rescaleCEF.py`: import-statement (line 1)
- `rescaleCEF.py`: function:rescaleCEF (line 23)
- `rescaleCEF.py`: function:rescaleCEF (line 23)
- `rescaleCEF.py`: function:rescaleCEF (line 23)
- `rescaleCEF.py`: function:rescaleCEF (line 23)
- `wybourne_stevens.py`: import-statement (line 1)
- `wybourne_stevens.py`: function:_convert_parameters (line 75)
</details>

#### **calculate_radial_integral_RE()** *(line 542)*

- **Used by:** 10 files

<details><summary>View all 10 usages</summary>

- `ligands.py`: import-statement (line 24)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 258)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 258)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 592)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 592)
- `rescaleCEF.py`: import-statement (line 1)
- `rescaleCEF.py`: function:rescaleCEF (line 23)
- `rescaleCEF.py`: function:rescaleCEF (line 23)
- `rescaleCEF.py`: function:rescaleCEF (line 23)
- `rescaleCEF.py`: function:rescaleCEF (line 23)
</details>

#### **calculate_radial_integral_TM()** *(line 556)*

- **Used by:** 3 files

<details><summary>View all 3 usages</summary>

- `ligands.py`: import-statement (line 24)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 687)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 687)
</details>

#### **is_half_filled()** *(line 602)*

- **Used by:** 3 files

<details><summary>View all 3 usages</summary>

- `ligands.py`: import-statement (line 24)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 638)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 638)
</details>

### 🔧 Constants

- **ION_NUMS_TRANS_METAL** *(line 19)* - 0 usages
- **ION_NUMS_RARE_EARTH** *(line 45)* - 11 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: import-statement (line 8)
  - `cf_levels.py`: class:CFLevels → method:Bdict (line 30)
  - `cf_levels.py`: function:LandeGFactor (line 754)
  - `form_factors.py`: import-statement (line 2)
  - `form_factors.py`: function:RE_FormFactor (line 27)
  - `import_CIF.py`: import-statement (line 6)
  - `import_CIF.py`: function:importCIF (line 83)
  - `ligands.py`: import-statement (line 24)
  - `ligands.py`: class:Ligands → method:PointChargeModel (line 233)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 467)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 468)
  </details>

- **WYBOURNE_STEVENS_CONSTS** *(line 59)* - 2 usages
  <details><summary>View usages</summary>

  - `wybourne_stevens.py`: import-statement (line 1)
  - `wybourne_stevens.py`: function:_convert_parameters (line 81)
  </details>

- **SPIN_ORBIT_COUPLING_CM** *(line 339)* - 0 usages
- **SPIN_ORBIT_COUPLING_CONSTANTS** *(line 418)* - 4 usages
  <details><summary>View usages</summary>

  - `import_CIF.py`: import-statement (line 6)
  - `import_CIF.py`: function:importCIF (line 106)
  - `import_CIF.py`: function:importCIF (line 110)
  - `import_CIF.py`: function:checkTMexist (line 134)
  </details>

- **RADIAL_INTEGRALS_RARE_EARTH** *(line 496)* - 0 usages
- **RADIAL_INTEGRALS_TRANS_METAL** *(line 513)* - 2 usages
  <details><summary>View usages</summary>

  - `import_CIF.py`: import-statement (line 6)
  - `import_CIF.py`: function:checkTMexist (line 132)
  </details>

- **ION_HALF_FILLED** *(line 568)* - 2 usages
  <details><summary>View usages</summary>

  - `import_CIF.py`: import-statement (line 6)
  - `import_CIF.py`: function:checkTMexist (line 136)
  </details>

- **ION_NOT_HALF_FILLED** *(line 581)* - 2 usages
  <details><summary>View usages</summary>

  - `import_CIF.py`: import-statement (line 6)
  - `import_CIF.py`: function:checkTMexist (line 136)
  </details>

### 📌 Variables

- **ahc** *(line 616)* - 8 usages
  <details><summary>View usages</summary>

  - `ligands.py`: class:Ligands → method:PointChargeModel (line 236)
  - `ligands.py`: class:Ligands → method:PointChargeModel (line 258)
  - `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 571)
  - `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 592)
  - `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 663)
  - `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 687)
  - `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 760)
  - `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 784)
  </details>

- **a0** *(line 617)* - 8 usages
  <details><summary>View usages</summary>

  - `ligands.py`: class:Ligands → method:PointChargeModel (line 237)
  - `ligands.py`: class:Ligands → method:PointChargeModel (line 258)
  - `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 572)
  - `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 592)
  - `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 664)
  - `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 687)
  - `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 761)
  - `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 784)
  </details>

- **muB** *(line 618)* - 21 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:magnetization (line 407)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 409)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 505)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 537)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 604)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 610)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 638)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetization (line 1057)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetization (line 1059)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1192)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1194)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1198)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1200)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1204)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1206)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1210)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1212)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1247)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1259)
  - `fundamental_operators.py`: class:LSOperator → method:magnetization (line 592)
  - `fundamental_operators.py`: class:LSOperator → method:magnetization (line 593)
  </details>

- **k_B** *(line 619)* - 32 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:magnetization (line 430)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 433)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 434)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 438)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 440)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 506)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 516)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 521)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 521)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 531)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 532)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 630)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 632)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 633)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetization (line 1078)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetization (line 1081)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetization (line 1082)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetization (line 1086)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetization (line 1088)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1241)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1245)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1246)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1251)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1253)
  - `fundamental_operators.py`: class:LSOperator → method:magnetization (line 611)
  - `fundamental_operators.py`: class:LSOperator → method:magnetization (line 614)
  - `fundamental_operators.py`: class:LSOperator → method:magnetization (line 615)
  - `fundamental_operators.py`: class:LSOperator → method:magnetization (line 619)
  - `fundamental_operators.py`: class:LSOperator → method:magnetization (line 621)
  - `thermo_functions.py`: import-statement (line 2)
  - `thermo_functions.py`: function:partition_func (line 18)
  - `thermo_functions.py`: function:Cp1T (line 34)
  </details>

---

## 📁 create_fit_function.py {#create_fit_function-py}

### 📊 Functions

#### **makeFitFunction()** *(line 7)*

- **Used by:** 12 files

<details><summary>View all 12 usages</summary>

- `cf_levels.py`: import-statement (line 10)
- `cf_levels.py`: class:CFLevels → method:fitdata (line 651)
- `cf_levels.py`: class:CFLevels → method:fitdata (line 651)
- `cf_levels.py`: class:CFLevels → method:fitdata_GlobalOpt (line 671)
- `cf_levels.py`: class:CFLevels → method:fitdata_GlobalOpt (line 671)
- `cf_levels.py`: class:LS_CFLevels → method:fitdata (line 1370)
- `cf_levels.py`: class:LS_CFLevels → method:fitdata (line 1370)
- `ligands.py`: import-statement (line 26)
- `ligands.py`: class:Ligands → method:FitCharges (line 322)
- `ligands.py`: class:Ligands → method:FitCharges (line 322)
- `ligands.py`: class:LS_Ligands → method:FitChargesNeutrons (line 827)
- `ligands.py`: class:LS_Ligands → method:FitChargesNeutrons (line 827)
</details>

#### **resultfunc()** *(line 51)*

- **Used by:** 0 file

#### **makeCurveFitFunction()** *(line 61)*

- **Used by:** 0 file

#### **resultfunc()** *(line 106)*

- **Used by:** 0 file

---

## 📁 form_factors.py {#form_factors-py}

### 📊 Functions

#### **importRE_FF()** *(line 4)*

- **Used by:** 0 file

#### **RE_FormFactor()** *(line 14)*

- **Used by:** 7 files

<details><summary>View all 7 usages</summary>

- `cf_levels.py`: import-statement (line 9)
- `cf_levels.py`: class:CFLevels → method:neutronSpectrum2D (line 271)
- `cf_levels.py`: class:CFLevels → method:neutronSpectrum2D (line 271)
- `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum2D (line 280)
- `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum2D (line 280)
- `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum2D (line 933)
- `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum2D (line 933)
</details>

---

## 📁 fundamental_operators.py {#fundamental_operators-py}

### 🏛️ Classes

#### **Ket** *(line 46)*

- **Used by:** 19 files


**Methods:**

- `__init__()` *(line 59)* - 0 usages
- `Jz()` *(line 71)* - 106 usages
  <details><summary>View all 106 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:__init__ (line 25)
  - `cf_levels.py`: class:CFLevels → method:Hamiltonian (line 49)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 295)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 295)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 381)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 400)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 409)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 429)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 524)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 524)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 536)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 536)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 572)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 574)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 575)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 576)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 577)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 600)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 600)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 610)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 629)
  - `cf_levels.py`: class:opttransition → method:__init__ (line 728)
  - `cf_levels.py`: class:opttransition → method:__init__ (line 731)
  - `cf_levels.py`: class:opttransition → method:transition (line 736)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 802)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1278)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1278)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1296)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1297)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1298)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1299)
  - `stevens_operators.py`: function:StevensOp (line 40)
  - `stevens_operators.py`: function:StevensOp (line 40)
  - `stevens_operators.py`: function:StevensOp (line 48)
  - `stevens_operators.py`: function:StevensOp (line 57)
  - `stevens_operators.py`: function:StevensOp (line 57)
  - `stevens_operators.py`: function:StevensOp (line 59)
  - `stevens_operators.py`: function:StevensOp (line 61)
  - `stevens_operators.py`: function:StevensOp (line 61)
  - `stevens_operators.py`: function:StevensOp (line 68)
  - `stevens_operators.py`: function:StevensOp (line 68)
  - `stevens_operators.py`: function:StevensOp (line 70)
  - `stevens_operators.py`: function:StevensOp (line 70)
  - `stevens_operators.py`: function:StevensOp (line 72)
  - `stevens_operators.py`: function:StevensOp (line 72)
  - `stevens_operators.py`: function:StevensOp (line 74)
  - `stevens_operators.py`: function:StevensOp (line 74)
  - `stevens_operators.py`: function:StevensOp (line 76)
  - `stevens_operators.py`: function:StevensOp (line 76)
  - `stevens_operators.py`: function:StevensOp (line 83)
  - `stevens_operators.py`: function:StevensOp (line 83)
  - `stevens_operators.py`: function:StevensOp (line 85)
  - `stevens_operators.py`: function:StevensOp (line 85)
  - `stevens_operators.py`: function:StevensOp (line 87)
  - `stevens_operators.py`: function:StevensOp (line 87)
  - `stevens_operators.py`: function:StevensOp (line 87)
  - `stevens_operators.py`: function:StevensOp (line 87)
  - `stevens_operators.py`: function:StevensOp (line 89)
  - `stevens_operators.py`: function:StevensOp (line 89)
  - `stevens_operators.py`: function:StevensOp (line 93)
  - `stevens_operators.py`: function:StevensOp (line 93)
  - `stevens_operators.py`: function:StevensOp (line 95)
  - `stevens_operators.py`: function:StevensOp (line 95)
  - `stevens_operators.py`: function:StevensOp (line 97)
  - `stevens_operators.py`: function:StevensOp (line 97)
  - `stevens_operators.py`: function:StevensOp (line 97)
  - `stevens_operators.py`: function:StevensOp (line 97)
  - `stevens_operators.py`: function:StevensOp (line 102)
  - `stevens_operators.py`: function:StevensOp (line 102)
  - `stevens_operators.py`: function:StevensOp (line 104)
  - `stevens_operators.py`: function:StevensOp (line 104)
  - `stevens_operators.py`: function:StevensOp (line 106)
  - `stevens_operators.py`: function:StevensOp (line 106)
  - `stevens_operators.py`: function:StevensOp (line 106)
  - `stevens_operators.py`: function:StevensOp (line 106)
  - `stevens_operators.py`: function:StevensOp (line 108)
  - `stevens_operators.py`: function:StevensOp (line 108)
  - `stevens_operators.py`: function:StevensOp (line 109)
  - `stevens_operators.py`: function:StevensOp (line 109)
  - `stevens_operators.py`: function:StevensOp (line 111)
  - `stevens_operators.py`: function:StevensOp (line 111)
  - `stevens_operators.py`: function:StevensOp (line 111)
  - `stevens_operators.py`: function:StevensOp (line 112)
  - `stevens_operators.py`: function:StevensOp (line 112)
  - `stevens_operators.py`: function:StevensOp (line 112)
  - `stevens_operators.py`: function:StevensOp (line 114)
  - `stevens_operators.py`: function:StevensOp (line 114)
  - `stevens_operators.py`: function:StevensOp (line 114)
  - `stevens_operators.py`: function:StevensOp (line 119)
  - `stevens_operators.py`: function:StevensOp (line 119)
  - `stevens_operators.py`: function:StevensOp (line 121)
  - `stevens_operators.py`: function:StevensOp (line 121)
  - `stevens_operators.py`: function:StevensOp (line 123)
  - `stevens_operators.py`: function:StevensOp (line 123)
  - `stevens_operators.py`: function:StevensOp (line 123)
  - `stevens_operators.py`: function:StevensOp (line 123)
  - `stevens_operators.py`: function:StevensOp (line 125)
  - `stevens_operators.py`: function:StevensOp (line 125)
  - `stevens_operators.py`: function:StevensOp (line 126)
  - `stevens_operators.py`: function:StevensOp (line 126)
  - `stevens_operators.py`: function:StevensOp (line 128)
  - `stevens_operators.py`: function:StevensOp (line 128)
  - `stevens_operators.py`: function:StevensOp (line 128)
  - `stevens_operators.py`: function:StevensOp (line 129)
  - `stevens_operators.py`: function:StevensOp (line 129)
  - `stevens_operators.py`: function:StevensOp (line 129)
  </details>

- `Jplus()` *(line 75)* - 5 usages
  <details><summary>View all 5 usages</summary>

  - `fundamental_operators.py`: class:Ket → method:Jx (line 99)
  - `fundamental_operators.py`: class:Ket → method:Jy (line 107)
  - `fundamental_operators.py`: class:Operator → method:Jx (line 267)
  - `fundamental_operators.py`: class:Operator → method:Jy (line 276)
  - `stevens_operators.py`: function:StevensOp (line 41)
  </details>

- `Jminus()` *(line 84)* - 5 usages
  <details><summary>View all 5 usages</summary>

  - `fundamental_operators.py`: class:Ket → method:Jx (line 99)
  - `fundamental_operators.py`: class:Ket → method:Jy (line 107)
  - `fundamental_operators.py`: class:Operator → method:Jx (line 268)
  - `fundamental_operators.py`: class:Operator → method:Jy (line 277)
  - `stevens_operators.py`: function:StevensOp (line 42)
  </details>

- `Jx()` *(line 93)* - 31 usages
  <details><summary>View all 31 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:__init__ (line 25)
  - `cf_levels.py`: class:CFLevels → method:Hamiltonian (line 49)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 293)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 293)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 379)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 398)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 409)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 427)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 522)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 522)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 534)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 534)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 570)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 580)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 581)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 582)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 583)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 598)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 598)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 610)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 627)
  - `cf_levels.py`: class:opttransition → method:__init__ (line 726)
  - `cf_levels.py`: class:opttransition → method:__init__ (line 729)
  - `cf_levels.py`: class:opttransition → method:transition (line 734)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 800)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1278)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1278)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1301)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1302)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1303)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1304)
  </details>

- `Jy()` *(line 101)* - 31 usages
  <details><summary>View all 31 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:__init__ (line 25)
  - `cf_levels.py`: class:CFLevels → method:Hamiltonian (line 49)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 294)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 294)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 380)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 399)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 409)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 428)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 523)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 523)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 535)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 535)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 571)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 585)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 586)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 587)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 588)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 599)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 599)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 610)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 628)
  - `cf_levels.py`: class:opttransition → method:__init__ (line 727)
  - `cf_levels.py`: class:opttransition → method:__init__ (line 730)
  - `cf_levels.py`: class:opttransition → method:transition (line 735)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 801)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1278)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1278)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1306)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1307)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1308)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1309)
  </details>

- `R()` *(line 109)* - 2 usages
  <details><summary>View usages</summary>

  - `thermo_functions.py`: function:Cp1T (line 33)
  - `thermo_functions.py`: function:Cp1T (line 38)
  </details>

- `_Rz()` *(line 123)* - 1 usages
  <details><summary>View usages</summary>

  - `fundamental_operators.py`: class:Ket → method:R (line 121)
  </details>

- `_Ry()` *(line 134)* - 0 usages
- `_WignersFormula()` *(line 148)* - 1 usages
  <details><summary>View usages</summary>

  - `fundamental_operators.py`: class:Ket → method:_Ry (line 145)
  </details>

- `__mul__()` *(line 178)* - 0 usages
- `__add__()` *(line 191)* - 0 usages

<details><summary>📍 View all 19 class usages</summary>

- `cf_levels.py`: import-statement (line 11)
- `cf_levels.py`: class:CFLevels → method:gsExpectation (line 378)
- `cf_levels.py`: class:CFLevels → method:gsExpectation (line 378)
- `cf_levels.py`: class:CFLevels → method:magnetization (line 424)
- `cf_levels.py`: class:CFLevels → method:magnetization (line 424)
- `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 520)
- `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 520)
- `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 530)
- `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 530)
- `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 626)
- `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 626)
- `cf_levels.py`: class:CFLevels → method:testEigenvectors (line 699)
- `cf_levels.py`: class:CFLevels → method:testEigenvectors (line 699)
- `cf_levels.py`: class:CFLevels → method:testEigenvectors (line 699)
- `cf_levels.py`: class:CFLevels → method:testEigenvectors (line 699)
- `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum (line 901)
- `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum (line 901)
- `cf_levels.py`: class:LS_CFLevels → method:normalizedNeutronSpectrum (line 950)
- `cf_levels.py`: class:LS_CFLevels → method:normalizedNeutronSpectrum (line 950)
</details>

#### **Operator** *(line 199)*

- **Used by:** 14 files


**Methods:**

- `__init__()` *(line 212)* - 0 usages
- `Jz()` *(line 221)* - 106 usages
  <details><summary>View all 106 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:__init__ (line 25)
  - `cf_levels.py`: class:CFLevels → method:Hamiltonian (line 49)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 295)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 295)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 381)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 400)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 409)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 429)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 524)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 524)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 536)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 536)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 572)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 574)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 575)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 576)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 577)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 600)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 600)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 610)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 629)
  - `cf_levels.py`: class:opttransition → method:__init__ (line 728)
  - `cf_levels.py`: class:opttransition → method:__init__ (line 731)
  - `cf_levels.py`: class:opttransition → method:transition (line 736)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 802)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1278)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1278)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1296)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1297)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1298)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1299)
  - `stevens_operators.py`: function:StevensOp (line 40)
  - `stevens_operators.py`: function:StevensOp (line 40)
  - `stevens_operators.py`: function:StevensOp (line 48)
  - `stevens_operators.py`: function:StevensOp (line 57)
  - `stevens_operators.py`: function:StevensOp (line 57)
  - `stevens_operators.py`: function:StevensOp (line 59)
  - `stevens_operators.py`: function:StevensOp (line 61)
  - `stevens_operators.py`: function:StevensOp (line 61)
  - `stevens_operators.py`: function:StevensOp (line 68)
  - `stevens_operators.py`: function:StevensOp (line 68)
  - `stevens_operators.py`: function:StevensOp (line 70)
  - `stevens_operators.py`: function:StevensOp (line 70)
  - `stevens_operators.py`: function:StevensOp (line 72)
  - `stevens_operators.py`: function:StevensOp (line 72)
  - `stevens_operators.py`: function:StevensOp (line 74)
  - `stevens_operators.py`: function:StevensOp (line 74)
  - `stevens_operators.py`: function:StevensOp (line 76)
  - `stevens_operators.py`: function:StevensOp (line 76)
  - `stevens_operators.py`: function:StevensOp (line 83)
  - `stevens_operators.py`: function:StevensOp (line 83)
  - `stevens_operators.py`: function:StevensOp (line 85)
  - `stevens_operators.py`: function:StevensOp (line 85)
  - `stevens_operators.py`: function:StevensOp (line 87)
  - `stevens_operators.py`: function:StevensOp (line 87)
  - `stevens_operators.py`: function:StevensOp (line 87)
  - `stevens_operators.py`: function:StevensOp (line 87)
  - `stevens_operators.py`: function:StevensOp (line 89)
  - `stevens_operators.py`: function:StevensOp (line 89)
  - `stevens_operators.py`: function:StevensOp (line 93)
  - `stevens_operators.py`: function:StevensOp (line 93)
  - `stevens_operators.py`: function:StevensOp (line 95)
  - `stevens_operators.py`: function:StevensOp (line 95)
  - `stevens_operators.py`: function:StevensOp (line 97)
  - `stevens_operators.py`: function:StevensOp (line 97)
  - `stevens_operators.py`: function:StevensOp (line 97)
  - `stevens_operators.py`: function:StevensOp (line 97)
  - `stevens_operators.py`: function:StevensOp (line 102)
  - `stevens_operators.py`: function:StevensOp (line 102)
  - `stevens_operators.py`: function:StevensOp (line 104)
  - `stevens_operators.py`: function:StevensOp (line 104)
  - `stevens_operators.py`: function:StevensOp (line 106)
  - `stevens_operators.py`: function:StevensOp (line 106)
  - `stevens_operators.py`: function:StevensOp (line 106)
  - `stevens_operators.py`: function:StevensOp (line 106)
  - `stevens_operators.py`: function:StevensOp (line 108)
  - `stevens_operators.py`: function:StevensOp (line 108)
  - `stevens_operators.py`: function:StevensOp (line 109)
  - `stevens_operators.py`: function:StevensOp (line 109)
  - `stevens_operators.py`: function:StevensOp (line 111)
  - `stevens_operators.py`: function:StevensOp (line 111)
  - `stevens_operators.py`: function:StevensOp (line 111)
  - `stevens_operators.py`: function:StevensOp (line 112)
  - `stevens_operators.py`: function:StevensOp (line 112)
  - `stevens_operators.py`: function:StevensOp (line 112)
  - `stevens_operators.py`: function:StevensOp (line 114)
  - `stevens_operators.py`: function:StevensOp (line 114)
  - `stevens_operators.py`: function:StevensOp (line 114)
  - `stevens_operators.py`: function:StevensOp (line 119)
  - `stevens_operators.py`: function:StevensOp (line 119)
  - `stevens_operators.py`: function:StevensOp (line 121)
  - `stevens_operators.py`: function:StevensOp (line 121)
  - `stevens_operators.py`: function:StevensOp (line 123)
  - `stevens_operators.py`: function:StevensOp (line 123)
  - `stevens_operators.py`: function:StevensOp (line 123)
  - `stevens_operators.py`: function:StevensOp (line 123)
  - `stevens_operators.py`: function:StevensOp (line 125)
  - `stevens_operators.py`: function:StevensOp (line 125)
  - `stevens_operators.py`: function:StevensOp (line 126)
  - `stevens_operators.py`: function:StevensOp (line 126)
  - `stevens_operators.py`: function:StevensOp (line 128)
  - `stevens_operators.py`: function:StevensOp (line 128)
  - `stevens_operators.py`: function:StevensOp (line 128)
  - `stevens_operators.py`: function:StevensOp (line 129)
  - `stevens_operators.py`: function:StevensOp (line 129)
  - `stevens_operators.py`: function:StevensOp (line 129)
  </details>

- `Jplus()` *(line 235)* - 5 usages
  <details><summary>View all 5 usages</summary>

  - `fundamental_operators.py`: class:Ket → method:Jx (line 99)
  - `fundamental_operators.py`: class:Ket → method:Jy (line 107)
  - `fundamental_operators.py`: class:Operator → method:Jx (line 267)
  - `fundamental_operators.py`: class:Operator → method:Jy (line 276)
  - `stevens_operators.py`: function:StevensOp (line 41)
  </details>

- `Jminus()` *(line 249)* - 5 usages
  <details><summary>View all 5 usages</summary>

  - `fundamental_operators.py`: class:Ket → method:Jx (line 99)
  - `fundamental_operators.py`: class:Ket → method:Jy (line 107)
  - `fundamental_operators.py`: class:Operator → method:Jx (line 268)
  - `fundamental_operators.py`: class:Operator → method:Jy (line 277)
  - `stevens_operators.py`: function:StevensOp (line 42)
  </details>

- `Jx()` *(line 263)* - 31 usages
  <details><summary>View all 31 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:__init__ (line 25)
  - `cf_levels.py`: class:CFLevels → method:Hamiltonian (line 49)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 293)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 293)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 379)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 398)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 409)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 427)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 522)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 522)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 534)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 534)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 570)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 580)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 581)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 582)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 583)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 598)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 598)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 610)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 627)
  - `cf_levels.py`: class:opttransition → method:__init__ (line 726)
  - `cf_levels.py`: class:opttransition → method:__init__ (line 729)
  - `cf_levels.py`: class:opttransition → method:transition (line 734)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 800)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1278)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1278)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1301)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1302)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1303)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1304)
  </details>

- `Jy()` *(line 272)* - 31 usages
  <details><summary>View all 31 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:__init__ (line 25)
  - `cf_levels.py`: class:CFLevels → method:Hamiltonian (line 49)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 294)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 294)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 380)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 399)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 409)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 428)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 523)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 523)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 535)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 535)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 571)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 585)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 586)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 587)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 588)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 599)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 599)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 610)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 628)
  - `cf_levels.py`: class:opttransition → method:__init__ (line 727)
  - `cf_levels.py`: class:opttransition → method:__init__ (line 730)
  - `cf_levels.py`: class:opttransition → method:transition (line 735)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 801)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1278)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1278)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1306)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1307)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1308)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1309)
  </details>

- `__add__()` *(line 280)* - 0 usages
- `__radd__()` *(line 289)* - 0 usages
- `__sub__()` *(line 298)* - 0 usages
- `__mul__()` *(line 307)* - 0 usages
- `__rmul__()` *(line 320)* - 0 usages
- `__pow__()` *(line 329)* - 0 usages
- `__neg__()` *(line 341)* - 0 usages
- `__repr__()` *(line 347)* - 0 usages

<details><summary>📍 View all 14 class usages</summary>

- `cf_levels.py`: import-statement (line 11)
- `cf_levels.py`: class:CFLevels → method:__init__ (line 25)
- `cf_levels.py`: class:CFLevels → method:__init__ (line 25)
- `cf_levels.py`: class:CFLevels → method:__init__ (line 25)
- `cf_levels.py`: class:CFLevels → method:Hamiltonian (line 49)
- `cf_levels.py`: class:CFLevels → method:Hamiltonian (line 49)
- `cf_levels.py`: class:CFLevels → method:Hamiltonian (line 49)
- `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 598)
- `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 599)
- `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 600)
- `stevens_operators.py`: import-statement (line 31)
- `stevens_operators.py`: function:StevensOp (line 40)
- `stevens_operators.py`: function:StevensOp (line 41)
- `stevens_operators.py`: function:StevensOp (line 42)
</details>

#### **LSOperator** *(line 352)*

- **Used by:** 25 files


**Methods:**

- `__init__()` *(line 369)* - 0 usages
- `Lz()` *(line 380)* - 17 usages
  <details><summary>View all 17 usages</summary>

  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 789)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 789)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 791)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 802)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 805)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1340)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1340)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1344)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1345)
  - `constants.py`: function:PFgamma (line 269)
  - `constants.py`: function:PFgamma (line 273)
  - `constants.py`: function:PFgamma (line 274)
  - `constants.py`: function:PFgamma (line 275)
  - `fundamental_operators.py`: class:LSOperator → method:magnetization (line 582)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 481)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 481)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 483)
  </details>

- `Lplus()` *(line 394)* - 2 usages
  <details><summary>View usages</summary>

  - `fundamental_operators.py`: class:LSOperator → method:Lx (line 424)
  - `fundamental_operators.py`: class:LSOperator → method:Ly (line 431)
  </details>

- `Lminus()` *(line 408)* - 2 usages
  <details><summary>View usages</summary>

  - `fundamental_operators.py`: class:LSOperator → method:Lx (line 425)
  - `fundamental_operators.py`: class:LSOperator → method:Ly (line 432)
  </details>

- `Lx()` *(line 422)* - 13 usages
  <details><summary>View all 13 usages</summary>

  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 787)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 787)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 791)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 800)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 803)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1338)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1338)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1344)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1345)
  - `fundamental_operators.py`: class:LSOperator → method:magnetization (line 580)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 479)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 479)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 483)
  </details>

- `Ly()` *(line 429)* - 13 usages
  <details><summary>View all 13 usages</summary>

  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 788)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 788)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 791)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 801)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 804)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1339)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1339)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1344)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1345)
  - `fundamental_operators.py`: class:LSOperator → method:magnetization (line 581)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 480)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 480)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 483)
  </details>

- `Sz()` *(line 439)* - 9 usages
  <details><summary>View all 9 usages</summary>

  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 786)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 786)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 791)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 802)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 805)
  - `fundamental_operators.py`: class:LSOperator → method:magnetization (line 585)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 478)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 478)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 483)
  </details>

- `Splus()` *(line 453)* - 2 usages
  <details><summary>View usages</summary>

  - `fundamental_operators.py`: class:LSOperator → method:Sx (line 483)
  - `fundamental_operators.py`: class:LSOperator → method:Sy (line 490)
  </details>

- `Sminus()` *(line 467)* - 2 usages
  <details><summary>View usages</summary>

  - `fundamental_operators.py`: class:LSOperator → method:Sx (line 484)
  - `fundamental_operators.py`: class:LSOperator → method:Sy (line 491)
  </details>

- `Sx()` *(line 481)* - 9 usages
  <details><summary>View all 9 usages</summary>

  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 784)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 784)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 791)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 800)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 803)
  - `fundamental_operators.py`: class:LSOperator → method:magnetization (line 583)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 476)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 476)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 483)
  </details>

- `Sy()` *(line 488)* - 9 usages
  <details><summary>View all 9 usages</summary>

  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 785)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 785)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 791)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 801)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 804)
  - `fundamental_operators.py`: class:LSOperator → method:magnetization (line 584)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 477)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 477)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 483)
  </details>

- `__add__()` *(line 495)* - 0 usages
- `__radd__()` *(line 504)* - 0 usages
- `__sub__()` *(line 513)* - 0 usages
- `__mul__()` *(line 522)* - 0 usages
- `__rmul__()` *(line 531)* - 0 usages
- `__pow__()` *(line 540)* - 0 usages
- `__neg__()` *(line 548)* - 0 usages
- `__repr__()` *(line 554)* - 0 usages
- `magnetization()` *(line 558)* - 48 usages
  <details><summary>View all 48 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 463)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 464)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 465)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 466)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 472)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 473)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 474)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 475)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 481)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 482)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 483)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 484)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 492)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 493)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 494)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 495)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1107)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1108)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1109)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1110)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1116)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1117)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1118)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1119)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1125)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1126)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1127)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1128)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1136)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1137)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1138)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1139)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 653)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 654)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 655)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 656)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 662)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 663)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 664)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 665)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 671)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 672)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 673)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 674)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 682)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 683)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 684)
  - `fundamental_operators.py`: class:LSOperator → method:susceptibility (line 685)
  </details>

- `susceptibility()` *(line 626)* - 0 usages

<details><summary>📍 View all 25 class usages</summary>

- `cf_levels.py`: import-statement (line 11)
- `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 773)
- `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 773)
- `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 784)
- `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 785)
- `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 786)
- `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 787)
- `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 788)
- `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 789)
- `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1338)
- `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1339)
- `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1340)
- `ligands.py`: import-statement (line 28)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 476)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 477)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 478)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 479)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 480)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 481)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 604)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 604)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 699)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 699)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 796)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 796)
</details>

---

## 📁 import_CIF.py {#import_CIF-py}

### 📊 Functions

#### **importCIF()** *(line 10)*

- **Used by:** 1 file

<details><summary>View all 1 usages</summary>

- `__init__.py`: import-statement (line 13)
</details>

#### **checkTMexist()** *(line 130)*

- **Used by:** 0 file

---

## 📁 lattice_class.py {#lattice_class-py}

### 🏛️ Classes

#### **lattice** *(line 15)*

- **Used by:** 12 files


**Methods:**

- `__init__()` *(line 35)* - 0 usages
- `reciplatt()` *(line 65)* - 1 usages
  <details><summary>View usages</summary>

  - `lattice_class.py`: class:lattice → method:__init__ (line 63)
  </details>

- `cartesian()` *(line 89)* - 0 usages
- `ABC()` *(line 117)* - 0 usages
- `inverseA()` *(line 138)* - 0 usages

<details><summary>📍 View all 12 class usages</summary>

- `cif_file.py`: import-statement (line 4)
- `cif_file.py`: class:CifFile → method:__init__ (line 152)
- `cif_file.py`: class:CifFile → method:__init__ (line 152)
- `ligands.py`: import-statement (line 23)
- `ligands.py`: class:Ligands → method:__init__ (line 74)
- `ligands.py`: class:Ligands → method:__init__ (line 74)
- `ligands.py`: class:Ligands → method:__init__ (line 78)
- `ligands.py`: class:Ligands → method:__init__ (line 78)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 455)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 455)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 459)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 459)
</details>

---

## 📁 ligands.py {#ligands-py}

### 🏛️ Classes

#### **Ligands** *(line 31)*

- **Used by:** 4 files


**Methods:**

- `__init__()` *(line 56)* - 0 usages
- `rotateLigands()` *(line 85)* - 0 usages
- `rotateLigandsZ()` *(line 104)* - 0 usages
- `_rotateMatrix()` *(line 118)* - 4 usages
  <details><summary>View all 4 usages</summary>

  - `ligands.py`: class:Ligands → method:rotateLigands (line 102)
  - `ligands.py`: class:Ligands → method:rotateLigandsZ (line 115)
  - `ligands.py`: class:LS_Ligands → method:rotateLigands (line 492)
  - `ligands.py`: class:LS_Ligands → method:rotateLigandsZ (line 497)
  </details>

- `exportCif()` *(line 158)* - 0 usages
- `PointChargeModel()` *(line 165)* - 5 usages
  <details><summary>View all 5 usages</summary>

  - `import_CIF.py`: function:importCIF (line 93)
  - `import_CIF.py`: function:importCIF (line 95)
  - `ligands.py`: class:Ligands → method:FitCharges (line 343)
  - `ligands.py`: class:Ligands → method:FitCharges (line 346)
  - `ligands.py`: class:LS_Ligands → method:FitChargesNeutrons (line 843)
  </details>

- `FitChargesNeutrons()` *(line 278)* - 0 usages
- `FitCharges()` *(line 286)* - 1 usages
  <details><summary>View usages</summary>

  - `ligands.py`: class:Ligands → method:FitChargesNeutrons (line 284)
  </details>


<details><summary>📍 View all 4 class usages</summary>

- `import_CIF.py`: import-statement (line 7)
- `import_CIF.py`: function:importCIF (line 88)
- `import_CIF.py`: function:importCIF (line 88)
- `__init__.py`: import-statement (line 11)
</details>

#### **LS_Ligands** *(line 404)*

- **Used by:** 8 files


**Methods:**

- `__init__()` *(line 436)* - 0 usages
- `rotateLigands()` *(line 488)* - 0 usages
- `rotateLigandsZ()` *(line 494)* - 0 usages
- `_rotateMatrix()` *(line 500)* - 4 usages
  <details><summary>View all 4 usages</summary>

  - `ligands.py`: class:Ligands → method:rotateLigands (line 102)
  - `ligands.py`: class:Ligands → method:rotateLigandsZ (line 115)
  - `ligands.py`: class:LS_Ligands → method:rotateLigands (line 492)
  - `ligands.py`: class:LS_Ligands → method:rotateLigandsZ (line 497)
  </details>

- `PointChargeModel()` *(line 524)* - 5 usages
  <details><summary>View all 5 usages</summary>

  - `import_CIF.py`: function:importCIF (line 93)
  - `import_CIF.py`: function:importCIF (line 95)
  - `ligands.py`: class:Ligands → method:FitCharges (line 343)
  - `ligands.py`: class:Ligands → method:FitCharges (line 346)
  - `ligands.py`: class:LS_Ligands → method:FitChargesNeutrons (line 843)
  </details>

- `TMPointChargeModel()` *(line 613)* - 2 usages
  <details><summary>View usages</summary>

  - `import_CIF.py`: function:importCIF (line 113)
  - `import_CIF.py`: function:importCIF (line 115)
  </details>

- `UnknownTMPointChargeModel()` *(line 710)* - 0 usages
- `FitChargesNeutrons()` *(line 805)* - 0 usages

<details><summary>📍 View all 8 class usages</summary>

- `import_CIF.py`: import-statement (line 7)
- `import_CIF.py`: function:importCIF (line 85)
- `import_CIF.py`: function:importCIF (line 85)
- `import_CIF.py`: function:importCIF (line 103)
- `import_CIF.py`: function:importCIF (line 103)
- `import_CIF.py`: function:importCIF (line 109)
- `import_CIF.py`: function:importCIF (line 109)
- `__init__.py`: import-statement (line 11)
</details>

### 📊 Functions

#### **exportLigandCif()** *(line 353)*

- **Used by:** 0 file

---

## 📁 moments_of_inertia.py {#moments_of_inertia-py}

### 📊 Functions

#### **MomIntertia()** *(line 5)*

- **Used by:** 0 file

#### **selectZaxisMI()** *(line 21)*

- **Used by:** 0 file

#### **ContinuousShapeMeasure()** *(line 40)*

- **Used by:** 0 file

#### **anglesToVector()** *(line 52)*

- **Used by:** 0 file

#### **rotationMatrix()** *(line 56)*

- **Used by:** 0 file

#### **rotateArbAxis()** *(line 72)*

- **Used by:** 0 file

#### **findZaxis_SOM_rotation()** *(line 80)*

- **Used by:** 0 file

#### **fitfun()** *(line 83)*

- **Used by:** 0 file

#### **findZaxis()** *(line 103)*

- **Used by:** 5 files

<details><summary>View all 5 usages</summary>

- `cif_symmetry_import.py`: import-statement (line 5)
- `cif_symmetry_import.py`: function:FindPointGroupSymOps (line 176)
- `cif_symmetry_import.py`: function:FindPointGroupSymOps (line 176)
- `cif_symmetry_import.py`: function:FindPointGroupSymOps (line 248)
- `cif_symmetry_import.py`: function:FindPointGroupSymOps (line 248)
</details>

---

## 📁 plotPCF.py {#plotPCF-py}

### 🏛️ Classes

#### **atomplot** *(line 43)*

- **Used by:** 0 file


**Methods:**

- `__init__()` *(line 63)* - 0 usages
- `plotatoms()` *(line 82)* - 1 usages
  <details><summary>View usages</summary>

  - `plotPCF.py`: function:plotPCF (line 32)
  </details>

- `plotaxes()` *(line 119)* - 1 usages
  <details><summary>View usages</summary>

  - `plotPCF.py`: function:plotPCF (line 33)
  </details>

- `plotabc()` *(line 139)* - 0 usages
- `_flatten()` *(line 161)* - 12 usages
  <details><summary>View all 12 usages</summary>

  - `plotPCF.py`: class:atomplot → method:plotaxes (line 130)
  - `plotPCF.py`: class:atomplot → method:plotaxes (line 131)
  - `plotPCF.py`: class:atomplot → method:plotaxes (line 132)
  - `plotPCF.py`: class:atomplot → method:plotaxes (line 135)
  - `plotPCF.py`: class:atomplot → method:plotaxes (line 136)
  - `plotPCF.py`: class:atomplot → method:plotaxes (line 137)
  - `plotPCF.py`: class:atomplot → method:plotabc (line 152)
  - `plotPCF.py`: class:atomplot → method:plotabc (line 153)
  - `plotPCF.py`: class:atomplot → method:plotabc (line 154)
  - `plotPCF.py`: class:atomplot → method:plotabc (line 157)
  - `plotPCF.py`: class:atomplot → method:plotabc (line 158)
  - `plotPCF.py`: class:atomplot → method:plotabc (line 159)
  </details>

### 📊 Functions

#### **plotPCF()** *(line 10)*

- **Used by:** 3 files

<details><summary>View all 3 usages</summary>

- `cif_symmetry_import.py`: import-statement (line 4)
- `cif_symmetry_import.py`: function:FindPointGroupSymOps (line 319)
- `cif_symmetry_import.py`: function:FindPointGroupSymOps (line 319)
</details>

---

## 📁 rescaleCEF.py {#rescaleCEF-py}

### 📊 Functions

#### **rescaleCEF()** *(line 5)*

- **Used by:** 0 file

---

## 📁 stevens_operators.py {#stevens_operators-py}

### 📊 Functions

#### **StevensOp()** *(line 35)*

- **Used by:** 14 files

<details><summary>View all 14 usages</summary>

- `cf_levels.py`: import-statement (line 12)
- `cf_levels.py`: class:CFLevels → method:Bdict (line 37)
- `cf_levels.py`: class:CFLevels → method:Bdict (line 37)
- `ligands.py`: import-statement (line 25)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 261)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 261)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 266)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 266)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 598)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 598)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 693)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 693)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 790)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 790)
</details>

#### **LS_StevensOp()** *(line 134)*

- **Used by:** 10 files

<details><summary>View all 10 usages</summary>

- `cf_levels.py`: import-statement (line 12)
- `cf_levels.py`: class:LS_CFLevels → method:Bdict (line 818)
- `cf_levels.py`: class:LS_CFLevels → method:Bdict (line 818)
- `ligands.py`: import-statement (line 25)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 595)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 595)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 690)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 690)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 787)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 787)
</details>

---

## 📁 thermo_functions.py {#thermo_functions-py}

### 📊 Functions

#### **partition_func()** *(line 6)*

- **Used by:** 0 file

#### **Cp_from_CEF()** *(line 20)*

- **Used by:** 0 file

#### **Cp1T()** *(line 32)*

- **Used by:** 0 file

---

## 📁 wybourne_stevens.py {#wybourne_stevens-py}

### 📊 Functions

#### **WybourneToStevens()** *(line 6)*

- **Used by:** 0 file

#### **StevensToWybourne()** *(line 34)*

- **Used by:** 0 file

#### **_convert_parameters()** *(line 60)*

- **Used by:** 0 file

---

## 📁 __init__.py {#__init__-py}

---

## 📈 Summary Statistics

- **Total files analyzed:** 17
- **Total entities discovered:** 208
- **Total cross-file dependencies:** 1064

### 🏆 Most Frequently Used Entities

| Entity | Type | File | Usage Count |
|--------|------|------|-------------|
| `Jz` | METHOD | fundamental_operators.py | 106 |
| `Jz` | METHOD | fundamental_operators.py | 106 |
| `magnetization` | METHOD | cf_levels.py | 48 |
| `magnetization` | METHOD | cf_levels.py | 48 |
| `magnetization` | METHOD | fundamental_operators.py | 48 |
| `k_B` | VARIABLE | constants.py | 32 |
| `Jx` | METHOD | fundamental_operators.py | 31 |
| `Jx` | METHOD | fundamental_operators.py | 31 |
| `Jy` | METHOD | fundamental_operators.py | 31 |
| `Jy` | METHOD | fundamental_operators.py | 31 |

### 🔗 Files with Most Dependencies

| File | Dependency Count |
|------|------------------|
| cf_levels.py | 522 |
| ligands.py | 161 |
| stevens_operators.py | 158 |
| fundamental_operators.py | 88 |
| import_CIF.py | 34 |

### ⚠️ Potentially Unused Entities (35 total)

<details><summary>View all unused entities</summary>

| Entity | Type | File | Line |
|--------|------|------|------|
| `spec` | VARIABLE | cf_levels.py | 717 |
| `opttransition` | CLASS | cf_levels.py | 724 |
| `LandeGFactor` | FUNCTION | cf_levels.py | 743 |
| `findRotationAxis` | FUNCTION | cif_symmetry_import.py | 332 |
| `makeSymOpMatrix` | FUNCTION | cif_symmetry_import.py | 382 |
| `ION_NUMS_TRANS_METAL` | CONSTANT | constants.py | 19 |
| `_tesseral_dispatch` | FUNCTION | constants.py | 117 |
| `PFgamma` | FUNCTION | constants.py | 238 |
| `SPIN_ORBIT_COUPLING_CM` | CONSTANT | constants.py | 339 |
| `RADIAL_INTEGRALS_RARE_EARTH` | CONSTANT | constants.py | 496 |
| `resultfunc` | FUNCTION | create_fit_function.py | 51 |
| `makeCurveFitFunction` | FUNCTION | create_fit_function.py | 61 |
| `resultfunc` | FUNCTION | create_fit_function.py | 106 |
| `importRE_FF` | FUNCTION | form_factors.py | 4 |
| `checkTMexist` | FUNCTION | import_CIF.py | 130 |
| `exportLigandCif` | FUNCTION | ligands.py | 353 |
| `MomIntertia` | FUNCTION | moments_of_inertia.py | 5 |
| `selectZaxisMI` | FUNCTION | moments_of_inertia.py | 21 |
| `ContinuousShapeMeasure` | FUNCTION | moments_of_inertia.py | 40 |
| `anglesToVector` | FUNCTION | moments_of_inertia.py | 52 |
| `rotationMatrix` | FUNCTION | moments_of_inertia.py | 56 |
| `rotateArbAxis` | FUNCTION | moments_of_inertia.py | 72 |
| `findZaxis_SOM_rotation` | FUNCTION | moments_of_inertia.py | 80 |
| `fitfun` | FUNCTION | moments_of_inertia.py | 83 |
| `atomplot` | CLASS | plotPCF.py | 43 |
| `rescaleCEF` | FUNCTION | rescaleCEF.py | 5 |
| `partition_func` | FUNCTION | thermo_functions.py | 6 |
| `Cp_from_CEF` | FUNCTION | thermo_functions.py | 20 |
| `Cp1T` | FUNCTION | thermo_functions.py | 32 |
| `WybourneToStevens` | FUNCTION | wybourne_stevens.py | 6 |
| `StevensToWybourne` | FUNCTION | wybourne_stevens.py | 34 |
| `_convert_parameters` | FUNCTION | wybourne_stevens.py | 60 |
| `__version__` | VARIABLE | __init__.py | 7 |
| `__author__` | VARIABLE | __init__.py | 8 |
| `__all__` | VARIABLE | __init__.py | 15 |

</details>