# 📊 Python Package Dependency Analysis Report

**Generated:** 2025-09-28 13:09:16

**Package:** `edited_package`

---

## 📑 Table of Contents

- [cf_levels.py](#cf_levels-py)
- [cif_file.py](#cif_file-py)
- [cif_symmetry_import.py](#cif_symmetry_import-py)
- [constants.py](#constants-py)
- [create_fit_function.py](#create_fit_function-py)
- [form_factors.py](#form_factors-py)
- [import_CIF.py](#import_CIF-py)
- [latex_cef_print.py](#latex_cef_print-py)
- [lattice_class.py](#lattice_class-py)
- [ligands.py](#ligands-py)
- [moments_of_inertia.py](#moments_of_inertia-py)
- [operators.py](#operators-py)
- [plot_ligands.py](#plot_ligands-py)
- [PyCrystalField.py](#PyCrystalField-py)
- [rescale_CEF.py](#rescale_CEF-py)
- [stevens_operators.py](#stevens_operators-py)
- [thermo_functions.py](#thermo_functions-py)
- [wybourne_stevens.py](#wybourne_stevens-py)
- [__init__.py](#__init__-py)
- [undefined_files\PCF_misc_functions.py](#undefined_files-PCF_misc_functions-py)

---

## 📁 cf_levels.py {#cf_levels-py}

### 🏛️ Classes

#### **CFLevels** *(line 14)*

- **Used by:** 2 files


**Methods:**

- `__init__()` *(line 17)* - 0 usages
- `Bdict()` *(line 39)* - 6 usages
  <details><summary>View all 6 usages</summary>

  - `wybourne_stevens.py`: function:WybourneToStevens (line 5)
  - `wybourne_stevens.py`: function:WybourneToStevens (line 9)
  - `wybourne_stevens.py`: function:WybourneToStevens (line 11)
  - `wybourne_stevens.py`: function:StevensToWybourne (line 16)
  - `wybourne_stevens.py`: function:StevensToWybourne (line 20)
  - `wybourne_stevens.py`: function:StevensToWybourne (line 22)
  </details>

- `Hamiltonian()` *(line 55)* - 4 usages
  <details><summary>View all 4 usages</summary>

  - `ligands.py`: class:Ligands → method:PointChargeModel (line 140)
  - `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 342)
  - `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 419)
  - `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 494)
  </details>

- `newCoeff()` *(line 62)* - 0 usages
- `diagonalize()` *(line 67)* - 14 usages
  <details><summary>View all 14 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:newCoeff (line 65)
  - `cf_levels.py`: class:CFLevels → method:neutronSpectrum (line 138)
  - `cf_levels.py`: class:CFLevels → method:neutronSpectrum_customLineshape (line 182)
  - `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum (line 213)
  - `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum_customLineshape (line 243)
  - `cf_levels.py`: class:CFLevels → method:printEigenvectors (line 329)
  - `cf_levels.py`: class:CFLevels → method:printLaTexEigenvectors (line 345)
  - `cf_levels.py`: class:LS_CFLevels → method:newCoeff (line 874)
  - `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum (line 913)
  - `cf_levels.py`: class:LS_CFLevels → method:normalizedNeutronSpectrum (line 960)
  - `cf_levels.py`: class:LS_CFLevels → method:printEigenvectors (line 1044)
  - `cf_levels.py`: class:LS_CFLevels → method:printLaTexEigenvectors (line 1403)
  - `ligands.py`: class:Ligands → method:FitCharges (line 182)
  - `ligands.py`: class:LS_Ligands → method:FitChargesNeutrons (line 521)
  </details>

- `diagonalize_banded()` *(line 88)* - 1 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:gtensor (line 538)
  </details>

- `_findbands()` *(line 106)* - 3 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:diagonalize (line 77)
  - `cf_levels.py`: class:CFLevels → method:diagonalize_banded (line 95)
  - `cf_levels.py`: class:LS_CFLevels → method:diagonalize (line 883)
  </details>

- `transitionIntensity()` *(line 119)* - 0 usages
- `neutronSpectrum()` *(line 130)* - 2 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:neutronSpectrum2D (line 265)
  - `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum2D (line 946)
  </details>

- `neutronSpectrum_customLineshape()` *(line 172)* - 0 usages
- `normalizedNeutronSpectrum()` *(line 205)* - 1 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum2D (line 274)
  </details>

- `normalizedNeutronSpectrum_customLineshape()` *(line 235)* - 0 usages
- `neutronSpectrum2D()` *(line 264)* - 0 usages
- `normalizedNeutronSpectrum2D()` *(line 273)* - 0 usages
- `_transition()` *(line 282)* - 3 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:testEigenvectors (line 691)
  - `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum (line 933)
  - `cf_levels.py`: class:LS_CFLevels → method:normalizedNeutronSpectrum (line 979)
  </details>

- `_lorentzian()` *(line 302)* - 0 usages
- `_voigt()` *(line 305)* - 4 usages
  <details><summary>View all 4 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:neutronSpectrum (line 156)
  - `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum (line 231)
  - `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum (line 936)
  - `cf_levels.py`: class:LS_CFLevels → method:normalizedNeutronSpectrum (line 982)
  </details>

- `_Re()` *(line 312)* - 14 usages
  <details><summary>View all 14 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:printEigenvectors (line 337)
  - `cf_levels.py`: class:CFLevels → method:printEigenvectors (line 337)
  - `cf_levels.py`: class:CFLevels → method:printLaTexEigenvectors (line 365)
  - `cf_levels.py`: class:CFLevels → method:printLaTexEigenvectors (line 366)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 377)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 378)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 379)
  - `cf_levels.py`: class:LS_CFLevels → method:printEigenvectors (line 1052)
  - `cf_levels.py`: class:LS_CFLevels → method:printEigenvectors (line 1052)
  - `cf_levels.py`: class:LS_CFLevels → method:gsExpectation (line 1061)
  - `cf_levels.py`: class:LS_CFLevels → method:gsExpectation (line 1062)
  - `cf_levels.py`: class:LS_CFLevels → method:gsExpectation (line 1063)
  - `cf_levels.py`: class:LS_CFLevels → method:printLaTexEigenvectors (line 1432)
  - `cf_levels.py`: class:LS_CFLevels → method:printLaTexEigenvectors (line 1433)
  </details>

- `printEigenvectors()` *(line 324)* - 0 usages
- `printLaTexEigenvectors()` *(line 340)* - 0 usages
- `gsExpectation()` *(line 370)* - 0 usages
- `magnetization()` *(line 383)* - 48 usages
  <details><summary>View all 48 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 459)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 460)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 461)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 462)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 468)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 469)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 470)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 471)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 477)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 478)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 479)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 480)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 488)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 489)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 490)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 491)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1122)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1123)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1124)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1125)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1131)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1132)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1133)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1134)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1140)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1141)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1142)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1143)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1151)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1152)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1153)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1154)
  - `operators.py`: class:LSOperator → method:susceptibility (line 423)
  - `operators.py`: class:LSOperator → method:susceptibility (line 424)
  - `operators.py`: class:LSOperator → method:susceptibility (line 425)
  - `operators.py`: class:LSOperator → method:susceptibility (line 426)
  - `operators.py`: class:LSOperator → method:susceptibility (line 432)
  - `operators.py`: class:LSOperator → method:susceptibility (line 433)
  - `operators.py`: class:LSOperator → method:susceptibility (line 434)
  - `operators.py`: class:LSOperator → method:susceptibility (line 435)
  - `operators.py`: class:LSOperator → method:susceptibility (line 441)
  - `operators.py`: class:LSOperator → method:susceptibility (line 442)
  - `operators.py`: class:LSOperator → method:susceptibility (line 443)
  - `operators.py`: class:LSOperator → method:susceptibility (line 444)
  - `operators.py`: class:LSOperator → method:susceptibility (line 452)
  - `operators.py`: class:LSOperator → method:susceptibility (line 453)
  - `operators.py`: class:LSOperator → method:susceptibility (line 454)
  - `operators.py`: class:LSOperator → method:susceptibility (line 455)
  </details>

- `susceptibility()` *(line 445)* - 0 usages
- `susceptibilityPert()` *(line 498)* - 0 usages
- `gtensor()` *(line 535)* - 0 usages
- `gtensorzeeman()` *(line 590)* - 0 usages
- `fitdata()` *(line 636)* - 0 usages
- `fitdata_GlobalOpt()` *(line 660)* - 0 usages
- `testEigenvectors()` *(line 680)* - 0 usages

<details><summary>📍 View all 2 class usages</summary>

- `ligands.py`: import-statement (line 9)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 140)
</details>

#### **OpticalTransition** *(line 710)*

- **Used by:** 0 file


**Methods:**

- `__init__()` *(line 721)* - 0 usages
- `transition_strength()` *(line 734)* - 0 usages
#### **LS_CFLevels** *(line 790)*

- **Used by:** 4 files


**Methods:**

- `__init__()` *(line 792)* - 0 usages
- `Bdict()` *(line 831)* - 6 usages
  <details><summary>View all 6 usages</summary>

  - `wybourne_stevens.py`: function:WybourneToStevens (line 5)
  - `wybourne_stevens.py`: function:WybourneToStevens (line 9)
  - `wybourne_stevens.py`: function:WybourneToStevens (line 11)
  - `wybourne_stevens.py`: function:StevensToWybourne (line 16)
  - `wybourne_stevens.py`: function:StevensToWybourne (line 20)
  - `wybourne_stevens.py`: function:StevensToWybourne (line 22)
  </details>

- `Hamiltonian()` *(line 846)* - 4 usages
  <details><summary>View all 4 usages</summary>

  - `ligands.py`: class:Ligands → method:PointChargeModel (line 140)
  - `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 342)
  - `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 419)
  - `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 494)
  </details>

- `newCoeff()` *(line 871)* - 0 usages
- `diagonalize()` *(line 876)* - 14 usages
  <details><summary>View all 14 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:newCoeff (line 65)
  - `cf_levels.py`: class:CFLevels → method:neutronSpectrum (line 138)
  - `cf_levels.py`: class:CFLevels → method:neutronSpectrum_customLineshape (line 182)
  - `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum (line 213)
  - `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum_customLineshape (line 243)
  - `cf_levels.py`: class:CFLevels → method:printEigenvectors (line 329)
  - `cf_levels.py`: class:CFLevels → method:printLaTexEigenvectors (line 345)
  - `cf_levels.py`: class:LS_CFLevels → method:newCoeff (line 874)
  - `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum (line 913)
  - `cf_levels.py`: class:LS_CFLevels → method:normalizedNeutronSpectrum (line 960)
  - `cf_levels.py`: class:LS_CFLevels → method:printEigenvectors (line 1044)
  - `cf_levels.py`: class:LS_CFLevels → method:printLaTexEigenvectors (line 1403)
  - `ligands.py`: class:Ligands → method:FitCharges (line 182)
  - `ligands.py`: class:LS_Ligands → method:FitChargesNeutrons (line 521)
  </details>

- `_findbands()` *(line 895)* - 3 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:diagonalize (line 77)
  - `cf_levels.py`: class:CFLevels → method:diagonalize_banded (line 95)
  - `cf_levels.py`: class:LS_CFLevels → method:diagonalize (line 883)
  </details>

- `neutronSpectrum()` *(line 908)* - 2 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:neutronSpectrum2D (line 265)
  - `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum2D (line 946)
  </details>

- `neutronSpectrum2D()` *(line 945)* - 0 usages
- `normalizedNeutronSpectrum()` *(line 954)* - 1 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum2D (line 274)
  </details>

- `_transition()` *(line 987)* - 3 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:testEigenvectors (line 691)
  - `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum (line 933)
  - `cf_levels.py`: class:LS_CFLevels → method:normalizedNeutronSpectrum (line 979)
  </details>

- `_lorentzian()` *(line 1014)* - 0 usages
- `_voigt()` *(line 1018)* - 4 usages
  <details><summary>View all 4 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:neutronSpectrum (line 156)
  - `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum (line 231)
  - `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum (line 936)
  - `cf_levels.py`: class:LS_CFLevels → method:normalizedNeutronSpectrum (line 982)
  </details>

- `_Re()` *(line 1026)* - 14 usages
  <details><summary>View all 14 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:printEigenvectors (line 337)
  - `cf_levels.py`: class:CFLevels → method:printEigenvectors (line 337)
  - `cf_levels.py`: class:CFLevels → method:printLaTexEigenvectors (line 365)
  - `cf_levels.py`: class:CFLevels → method:printLaTexEigenvectors (line 366)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 377)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 378)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 379)
  - `cf_levels.py`: class:LS_CFLevels → method:printEigenvectors (line 1052)
  - `cf_levels.py`: class:LS_CFLevels → method:printEigenvectors (line 1052)
  - `cf_levels.py`: class:LS_CFLevels → method:gsExpectation (line 1061)
  - `cf_levels.py`: class:LS_CFLevels → method:gsExpectation (line 1062)
  - `cf_levels.py`: class:LS_CFLevels → method:gsExpectation (line 1063)
  - `cf_levels.py`: class:LS_CFLevels → method:printLaTexEigenvectors (line 1432)
  - `cf_levels.py`: class:LS_CFLevels → method:printLaTexEigenvectors (line 1433)
  </details>

- `printEigenvectors()` *(line 1039)* - 0 usages
- `gsExpectation()` *(line 1055)* - 0 usages
- `magnetization()` *(line 1067)* - 48 usages
  <details><summary>View all 48 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 459)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 460)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 461)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 462)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 468)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 469)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 470)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 471)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 477)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 478)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 479)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 480)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 488)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 489)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 490)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 491)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1122)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1123)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1124)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1125)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1131)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1132)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1133)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1134)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1140)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1141)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1142)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1143)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1151)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1152)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1153)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1154)
  - `operators.py`: class:LSOperator → method:susceptibility (line 423)
  - `operators.py`: class:LSOperator → method:susceptibility (line 424)
  - `operators.py`: class:LSOperator → method:susceptibility (line 425)
  - `operators.py`: class:LSOperator → method:susceptibility (line 426)
  - `operators.py`: class:LSOperator → method:susceptibility (line 432)
  - `operators.py`: class:LSOperator → method:susceptibility (line 433)
  - `operators.py`: class:LSOperator → method:susceptibility (line 434)
  - `operators.py`: class:LSOperator → method:susceptibility (line 435)
  - `operators.py`: class:LSOperator → method:susceptibility (line 441)
  - `operators.py`: class:LSOperator → method:susceptibility (line 442)
  - `operators.py`: class:LSOperator → method:susceptibility (line 443)
  - `operators.py`: class:LSOperator → method:susceptibility (line 444)
  - `operators.py`: class:LSOperator → method:susceptibility (line 452)
  - `operators.py`: class:LSOperator → method:susceptibility (line 453)
  - `operators.py`: class:LSOperator → method:susceptibility (line 454)
  - `operators.py`: class:LSOperator → method:susceptibility (line 455)
  </details>

- `susceptibility()` *(line 1112)* - 0 usages
- `susceptibilityDeriv()` *(line 1161)* - 0 usages
- `magnetizationDeriv()` *(line 1200)* - 12 usages
  <details><summary>View all 12 usages</summary>

  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1171)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1172)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1173)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1174)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1180)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1181)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1182)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1183)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1189)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1190)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1191)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibilityDeriv (line 1192)
  </details>

- `gtensor()` *(line 1276)* - 0 usages
- `gtensorperturb()` *(line 1340)* - 0 usages
- `fitdata()` *(line 1370)* - 0 usages
- `printLaTexEigenvectors()` *(line 1398)* - 0 usages

<details><summary>📍 View all 4 class usages</summary>

- `ligands.py`: import-statement (line 9)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 342)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 419)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 494)
</details>

### 📊 Functions

#### **LandeGFactor()** *(line 765)*

- **Used by:** 0 file

---

## 📁 cif_file.py {#cif_file-py}

### 🏛️ Classes

#### **CifFile** *(line 9)*

- **Used by:** 3 files


**Methods:**

- `__init__()` *(line 11)* - 0 usages
- `SymOperate()` *(line 158)* - 2 usages
  <details><summary>View usages</summary>

  - `cif_file.py`: class:CifFile → method:MakeUnitCell (line 185)
  - `cif_symmetry_import.py`: function:FindPointGroupSymOps (line 29)
  </details>

- `MakeUnitCell()` *(line 180)* - 1 usages
  <details><summary>View usages</summary>

  - `cif_file.py`: class:CifFile → method:__init__ (line 150)
  </details>

- `StructureFactor()` *(line 218)* - 0 usages
- `MultipleScattering()` *(line 260)* - 0 usages
- `_destringify()` *(line 305)* - 11 usages
  <details><summary>View all 11 usages</summary>

  - `cif_file.py`: class:CifFile → method:__init__ (line 32)
  - `cif_file.py`: class:CifFile → method:__init__ (line 34)
  - `cif_file.py`: class:CifFile → method:__init__ (line 36)
  - `cif_file.py`: class:CifFile → method:__init__ (line 38)
  - `cif_file.py`: class:CifFile → method:__init__ (line 40)
  - `cif_file.py`: class:CifFile → method:__init__ (line 42)
  - `cif_file.py`: class:CifFile → method:__init__ (line 86)
  - `cif_file.py`: class:CifFile → method:__init__ (line 87)
  - `cif_file.py`: class:CifFile → method:__init__ (line 88)
  - `cif_file.py`: class:CifFile → method:__init__ (line 90)
  - `cif_file.py`: class:CifFile → method:__init__ (line 92)
  </details>

- `_defractionify()` *(line 323)* - 3 usages
  <details><summary>View usages</summary>

  - `cif_file.py`: class:CifFile → method:SymOperate (line 167)
  - `cif_file.py`: class:CifFile → method:SymOperate (line 168)
  - `cif_file.py`: class:CifFile → method:SymOperate (line 169)
  </details>

- `_duplicaterow()` *(line 327)* - 1 usages
  <details><summary>View usages</summary>

  - `cif_file.py`: class:CifFile → method:MakeUnitCell (line 193)
  </details>

- `_NumElements()` *(line 340)* - 1 usages
  <details><summary>View usages</summary>

  - `cif_file.py`: class:CifFile → method:StructureFactor (line 224)
  </details>

- `_kvector()` *(line 347)* - 1 usages
  <details><summary>View usages</summary>

  - `cif_file.py`: class:CifFile → method:MultipleScattering (line 263)
  </details>


<details><summary>📍 View all 3 class usages</summary>

- `import_CIF.py`: import-statement (line 5)
- `import_CIF.py`: function:importCIF (line 18)
- `import_CIF.py`: function:importCIF (line 18)
</details>

---

## 📁 cif_symmetry_import.py {#cif_symmetry_import-py}

### 📊 Functions

#### **FindPointGroupSymOps()** *(line 6)*

- **Used by:** 3 files

<details><summary>View all 3 usages</summary>

- `import_CIF.py`: import-statement (line 4)
- `import_CIF.py`: function:importCIF (line 66)
- `import_CIF.py`: function:importCIF (line 66)
</details>

#### **findRotationAxis()** *(line 332)*

- **Used by:** 0 file

#### **makeSymOpMatrix()** *(line 376)*

- **Used by:** 0 file

---

## 📁 constants.py {#constants-py}

### 📊 Functions

#### **Constant()** *(line 100)*

- **Used by:** 9 files

<details><summary>View all 9 usages</summary>

- `ligands.py`: import-statement (line 6)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 128)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 128)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 327)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 327)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 403)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 403)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 479)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 479)
</details>

#### **calculate_tesseral_harmonic()** *(line 105)*

- **Used by:** 9 files

<details><summary>View all 9 usages</summary>

- `ligands.py`: import-statement (line 6)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 124)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 124)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 323)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 323)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 399)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 399)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 475)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 475)
</details>

#### **_tesseral_dispatch()** *(line 117)*

- **Used by:** 0 file

#### **PFalpha()** *(line 189)*

- **Used by:** 5 files

<details><summary>View all 5 usages</summary>

- `ligands.py`: import-statement (line 6)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 387)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 387)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 463)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 463)
</details>

#### **PFbeta()** *(line 209)*

- **Used by:** 5 files

<details><summary>View all 5 usages</summary>

- `ligands.py`: import-statement (line 6)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 388)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 388)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 464)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 464)
</details>

#### **PFgamma()** *(line 240)*

- **Used by:** 0 file

#### **LStheta()** *(line 290)*

- **Used by:** 8 files

<details><summary>View all 8 usages</summary>

- `ligands.py`: import-statement (line 6)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 327)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 327)
- `wybourne_stevens.py`: import-statement (line 1)
- `wybourne_stevens.py`: function:WybourneToStevens (line 9)
- `wybourne_stevens.py`: function:WybourneToStevens (line 9)
- `wybourne_stevens.py`: function:StevensToWybourne (line 20)
- `wybourne_stevens.py`: function:StevensToWybourne (line 20)
</details>

#### **theta()** *(line 312)*

- **Used by:** 26 files

<details><summary>View all 26 usages</summary>

- `ligands.py`: import-statement (line 6)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 128)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 128)
- `moments_of_inertia.py`: function:anglesToVector (line 53)
- `moments_of_inertia.py`: function:anglesToVector (line 53)
- `moments_of_inertia.py`: function:anglesToVector (line 53)
- `moments_of_inertia.py`: function:rotationMatrix (line 58)
- `moments_of_inertia.py`: function:rotateArbAxis (line 73)
- `moments_of_inertia.py`: function:fitfun (line 84)
- `moments_of_inertia.py`: function:fitfun (line 85)
- `operators.py`: class:Ket → method:_Rz (line 34)
- `plot_ligands.py`: class:atomplot → method:__init__ (line 55)
- `plot_ligands.py`: class:atomplot → method:__init__ (line 55)
- `plot_ligands.py`: class:atomplot → method:__init__ (line 58)
- `plot_ligands.py`: class:atomplot → method:__init__ (line 58)
- `plot_ligands.py`: class:atomplot → method:__init__ (line 58)
- `rescale_CEF.py`: import-statement (line 1)
- `rescale_CEF.py`: function:rescaleCEF (line 7)
- `rescale_CEF.py`: function:rescaleCEF (line 7)
- `rescale_CEF.py`: function:rescaleCEF (line 7)
- `rescale_CEF.py`: function:rescaleCEF (line 7)
- `wybourne_stevens.py`: import-statement (line 1)
- `wybourne_stevens.py`: function:WybourneToStevens (line 11)
- `wybourne_stevens.py`: function:WybourneToStevens (line 11)
- `wybourne_stevens.py`: function:StevensToWybourne (line 22)
- `wybourne_stevens.py`: function:StevensToWybourne (line 22)
</details>

#### **calculate_radial_integral_RE()** *(line 544)*

- **Used by:** 10 files

<details><summary>View all 10 usages</summary>

- `ligands.py`: import-statement (line 6)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 128)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 128)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 327)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 327)
- `rescale_CEF.py`: import-statement (line 1)
- `rescale_CEF.py`: function:rescaleCEF (line 7)
- `rescale_CEF.py`: function:rescaleCEF (line 7)
- `rescale_CEF.py`: function:rescaleCEF (line 7)
- `rescale_CEF.py`: function:rescaleCEF (line 7)
</details>

#### **calculate_radial_integral_TM()** *(line 560)*

- **Used by:** 3 files

<details><summary>View all 3 usages</summary>

- `ligands.py`: import-statement (line 6)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 403)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 403)
</details>

#### **is_half_filled()** *(line 608)*

- **Used by:** 3 files

<details><summary>View all 3 usages</summary>

- `ligands.py`: import-statement (line 6)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 354)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 354)
</details>

### 🔧 Constants

- **ION_NUMS_TRANS_METAL** *(line 19)* - 0 usages
- **ION_NUMS_RARE_EARTH** *(line 45)* - 11 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: import-statement (line 7)
  - `cf_levels.py`: class:CFLevels → method:Bdict (line 40)
  - `cf_levels.py`: function:LandeGFactor (line 776)
  - `form_factors.py`: import-statement (line 2)
  - `form_factors.py`: function:RE_FormFactor (line 27)
  - `import_CIF.py`: import-statement (line 6)
  - `import_CIF.py`: function:importCIF (line 71)
  - `ligands.py`: import-statement (line 6)
  - `ligands.py`: class:Ligands → method:PointChargeModel (line 103)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 219)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 220)
  </details>

- **WYBOURNE_STEVENS_CONSTS** *(line 59)* - 5 usages
  <details><summary>View usages</summary>

  - `wybourne_stevens.py`: import-statement (line 1)
  - `wybourne_stevens.py`: function:WybourneToStevens (line 9)
  - `wybourne_stevens.py`: function:WybourneToStevens (line 11)
  - `wybourne_stevens.py`: function:StevensToWybourne (line 20)
  - `wybourne_stevens.py`: function:StevensToWybourne (line 22)
  </details>

- **SPIN_ORBIT_COUPLING_CM** *(line 341)* - 0 usages
- **SPIN_ORBIT_COUPLING_CONSTANTS** *(line 420)* - 4 usages
  <details><summary>View usages</summary>

  - `import_CIF.py`: import-statement (line 6)
  - `import_CIF.py`: function:importCIF (line 94)
  - `import_CIF.py`: function:importCIF (line 98)
  - `import_CIF.py`: function:checkTMexist (line 121)
  </details>

- **RADIAL_INTEGRALS_RARE_EARTH** *(line 498)* - 1 usages
  <details><summary>View usages</summary>

  - `constants.py`: import-statement (line 548)
  </details>

- **RADIAL_INTEGRALS_TRANS_METAL** *(line 515)* - 3 usages
  <details><summary>View usages</summary>

  - `constants.py`: import-statement (line 562)
  - `import_CIF.py`: import-statement (line 6)
  - `import_CIF.py`: function:checkTMexist (line 119)
  </details>

- **ION_HALF_FILLED** *(line 574)* - 3 usages
  <details><summary>View usages</summary>

  - `constants.py`: import-statement (line 611)
  - `import_CIF.py`: import-statement (line 6)
  - `import_CIF.py`: function:checkTMexist (line 123)
  </details>

- **ION_NOT_HALF_FILLED** *(line 587)* - 3 usages
  <details><summary>View usages</summary>

  - `constants.py`: import-statement (line 611)
  - `import_CIF.py`: import-statement (line 6)
  - `import_CIF.py`: function:checkTMexist (line 123)
  </details>

### 📌 Variables

- **ahc** *(line 624)* - 8 usages
  <details><summary>View usages</summary>

  - `ligands.py`: class:Ligands → method:PointChargeModel (line 106)
  - `ligands.py`: class:Ligands → method:PointChargeModel (line 128)
  - `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 306)
  - `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 327)
  - `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 379)
  - `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 403)
  - `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 455)
  - `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 479)
  </details>

- **a0** *(line 625)* - 8 usages
  <details><summary>View usages</summary>

  - `ligands.py`: class:Ligands → method:PointChargeModel (line 107)
  - `ligands.py`: class:Ligands → method:PointChargeModel (line 128)
  - `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 307)
  - `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 327)
  - `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 380)
  - `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 403)
  - `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 456)
  - `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 479)
  </details>

- **muB** *(line 626)* - 21 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:magnetization (line 404)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 406)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 501)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 533)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 598)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 604)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 632)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetization (line 1073)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetization (line 1075)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1207)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1209)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1213)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1215)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1219)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1221)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1225)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1227)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1262)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1274)
  - `operators.py`: class:LSOperator → method:magnetization (line 373)
  - `operators.py`: class:LSOperator → method:magnetization (line 374)
  </details>

- **k_B** *(line 627)* - 29 usages
  <details><summary>View usages</summary>

  - `cf_levels.py`: class:CFLevels → method:magnetization (line 427)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 430)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 431)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 435)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 437)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 502)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 512)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 517)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 517)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 527)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 528)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 624)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 626)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 627)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetization (line 1094)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetization (line 1097)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetization (line 1098)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetization (line 1102)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetization (line 1104)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1256)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1260)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1261)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1266)
  - `cf_levels.py`: class:LS_CFLevels → method:magnetizationDeriv (line 1268)
  - `operators.py`: class:LSOperator → method:magnetization (line 394)
  - `operators.py`: class:LSOperator → method:magnetization (line 397)
  - `operators.py`: class:LSOperator → method:magnetization (line 398)
  - `operators.py`: class:LSOperator → method:magnetization (line 402)
  - `operators.py`: class:LSOperator → method:magnetization (line 404)
  </details>

---

## 📁 create_fit_function.py {#create_fit_function-py}

### 📊 Functions

#### **makeFitFunction()** *(line 7)*

- **Used by:** 12 files

<details><summary>View all 12 usages</summary>

- `cf_levels.py`: import-statement (line 9)
- `cf_levels.py`: class:CFLevels → method:fitdata (line 644)
- `cf_levels.py`: class:CFLevels → method:fitdata (line 644)
- `cf_levels.py`: class:CFLevels → method:fitdata_GlobalOpt (line 664)
- `cf_levels.py`: class:CFLevels → method:fitdata_GlobalOpt (line 664)
- `cf_levels.py`: class:LS_CFLevels → method:fitdata (line 1381)
- `cf_levels.py`: class:LS_CFLevels → method:fitdata (line 1381)
- `ligands.py`: import-statement (line 8)
- `ligands.py`: class:Ligands → method:FitCharges (line 157)
- `ligands.py`: class:Ligands → method:FitCharges (line 157)
- `ligands.py`: class:LS_Ligands → method:FitChargesNeutrons (line 504)
- `ligands.py`: class:LS_Ligands → method:FitChargesNeutrons (line 504)
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

- `cf_levels.py`: import-statement (line 8)
- `cf_levels.py`: class:CFLevels → method:neutronSpectrum2D (line 270)
- `cf_levels.py`: class:CFLevels → method:neutronSpectrum2D (line 270)
- `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum2D (line 279)
- `cf_levels.py`: class:CFLevels → method:normalizedNeutronSpectrum2D (line 279)
- `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum2D (line 951)
- `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum2D (line 951)
</details>

---

## 📁 import_CIF.py {#import_CIF-py}

### 📊 Functions

#### **importCIF()** *(line 12)*

- **Used by:** 0 file

#### **checkTMexist()** *(line 118)*

- **Used by:** 0 file

---

## 📁 latex_cef_print.py {#latex_cef_print-py}

### 📊 Functions

#### **printLaTexCEFparams()** *(line 3)*

- **Used by:** 0 file

---

## 📁 lattice_class.py {#lattice_class-py}

### 🏛️ Classes

#### **lattice** *(line 3)*

- **Used by:** 12 files


**Methods:**

- `__init__()` *(line 4)* - 0 usages
- `reciplatt()` *(line 24)* - 1 usages
  <details><summary>View usages</summary>

  - `lattice_class.py`: class:lattice → method:__init__ (line 22)
  </details>

- `cartesian()` *(line 39)* - 0 usages
- `ABC()` *(line 53)* - 0 usages
- `inverseA()` *(line 62)* - 0 usages

<details><summary>📍 View all 12 class usages</summary>

- `cif_file.py`: import-statement (line 6)
- `cif_file.py`: class:CifFile → method:__init__ (line 154)
- `cif_file.py`: class:CifFile → method:__init__ (line 154)
- `ligands.py`: import-statement (line 4)
- `ligands.py`: class:Ligands → method:__init__ (line 19)
- `ligands.py`: class:Ligands → method:__init__ (line 19)
- `ligands.py`: class:Ligands → method:__init__ (line 23)
- `ligands.py`: class:Ligands → method:__init__ (line 23)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 207)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 207)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 211)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 211)
</details>

---

## 📁 ligands.py {#ligands-py}

### 🏛️ Classes

#### **Ligands** *(line 13)*

- **Used by:** 3 files


**Methods:**

- `__init__()` *(line 15)* - 0 usages
- `rotateLigands()` *(line 30)* - 0 usages
- `rotateLigandsZ()` *(line 36)* - 0 usages
- `_rotateMatrix()` *(line 43)* - 4 usages
  <details><summary>View all 4 usages</summary>

  - `ligands.py`: class:Ligands → method:rotateLigands (line 34)
  - `ligands.py`: class:Ligands → method:rotateLigandsZ (line 40)
  - `ligands.py`: class:LS_Ligands → method:rotateLigands (line 244)
  - `ligands.py`: class:LS_Ligands → method:rotateLigandsZ (line 250)
  </details>

- `exportCif()` *(line 67)* - 0 usages
- `PointChargeModel()` *(line 71)* - 5 usages
  <details><summary>View all 5 usages</summary>

  - `import_CIF.py`: function:importCIF (line 81)
  - `import_CIF.py`: function:importCIF (line 83)
  - `ligands.py`: class:Ligands → method:FitCharges (line 178)
  - `ligands.py`: class:Ligands → method:FitCharges (line 181)
  - `ligands.py`: class:LS_Ligands → method:FitChargesNeutrons (line 520)
  </details>

- `FitChargesNeutrons()` *(line 148)* - 0 usages
- `FitCharges()` *(line 153)* - 1 usages
  <details><summary>View usages</summary>

  - `ligands.py`: class:Ligands → method:FitChargesNeutrons (line 151)
  </details>


<details><summary>📍 View all 3 class usages</summary>

- `import_CIF.py`: import-statement (line 7)
- `import_CIF.py`: function:importCIF (line 77)
- `import_CIF.py`: function:importCIF (line 77)
</details>

#### **LS_Ligands** *(line 199)*

- **Used by:** 7 files


**Methods:**

- `__init__()` *(line 201)* - 0 usages
- `rotateLigands()` *(line 240)* - 0 usages
- `rotateLigandsZ()` *(line 246)* - 0 usages
- `_rotateMatrix()` *(line 253)* - 4 usages
  <details><summary>View all 4 usages</summary>

  - `ligands.py`: class:Ligands → method:rotateLigands (line 34)
  - `ligands.py`: class:Ligands → method:rotateLigandsZ (line 40)
  - `ligands.py`: class:LS_Ligands → method:rotateLigands (line 244)
  - `ligands.py`: class:LS_Ligands → method:rotateLigandsZ (line 250)
  </details>

- `PointChargeModel()` *(line 278)* - 5 usages
  <details><summary>View all 5 usages</summary>

  - `import_CIF.py`: function:importCIF (line 81)
  - `import_CIF.py`: function:importCIF (line 83)
  - `ligands.py`: class:Ligands → method:FitCharges (line 178)
  - `ligands.py`: class:Ligands → method:FitCharges (line 181)
  - `ligands.py`: class:LS_Ligands → method:FitChargesNeutrons (line 520)
  </details>

- `TMPointChargeModel()` *(line 348)* - 2 usages
  <details><summary>View usages</summary>

  - `import_CIF.py`: function:importCIF (line 101)
  - `import_CIF.py`: function:importCIF (line 103)
  </details>

- `UnknownTMPointChargeModel()` *(line 426)* - 0 usages
- `FitChargesNeutrons()` *(line 500)* - 0 usages

<details><summary>📍 View all 7 class usages</summary>

- `import_CIF.py`: import-statement (line 7)
- `import_CIF.py`: function:importCIF (line 73)
- `import_CIF.py`: function:importCIF (line 73)
- `import_CIF.py`: function:importCIF (line 91)
- `import_CIF.py`: function:importCIF (line 91)
- `import_CIF.py`: function:importCIF (line 97)
- `import_CIF.py`: function:importCIF (line 97)
</details>

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

- `cif_symmetry_import.py`: import-statement (line 4)
- `cif_symmetry_import.py`: function:FindPointGroupSymOps (line 167)
- `cif_symmetry_import.py`: function:FindPointGroupSymOps (line 167)
- `cif_symmetry_import.py`: function:FindPointGroupSymOps (line 241)
- `cif_symmetry_import.py`: function:FindPointGroupSymOps (line 241)
</details>

---

## 📁 operators.py {#operators-py}

### 🏛️ Classes

#### **Ket** *(line 3)*

- **Used by:** 19 files


**Methods:**

- `__init__()` *(line 4)* - 0 usages
- `Jz()` *(line 11)* - 105 usages
  <details><summary>View all 105 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:__init__ (line 32)
  - `cf_levels.py`: class:CFLevels → method:Hamiltonian (line 59)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 293)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 293)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 379)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 397)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 406)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 426)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 520)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 520)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 532)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 532)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 567)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 569)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 570)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 571)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 572)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 594)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 594)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 604)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 623)
  - `cf_levels.py`: class:OpticalTransition → method:__init__ (line 732)
  - `cf_levels.py`: class:OpticalTransition → method:transition_strength (line 751)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 824)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1290)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1290)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1308)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1309)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1310)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1311)
  - `stevens_operators.py`: function:StevensOp (line 8)
  - `stevens_operators.py`: function:StevensOp (line 8)
  - `stevens_operators.py`: function:StevensOp (line 16)
  - `stevens_operators.py`: function:StevensOp (line 25)
  - `stevens_operators.py`: function:StevensOp (line 25)
  - `stevens_operators.py`: function:StevensOp (line 27)
  - `stevens_operators.py`: function:StevensOp (line 29)
  - `stevens_operators.py`: function:StevensOp (line 29)
  - `stevens_operators.py`: function:StevensOp (line 36)
  - `stevens_operators.py`: function:StevensOp (line 36)
  - `stevens_operators.py`: function:StevensOp (line 38)
  - `stevens_operators.py`: function:StevensOp (line 38)
  - `stevens_operators.py`: function:StevensOp (line 40)
  - `stevens_operators.py`: function:StevensOp (line 40)
  - `stevens_operators.py`: function:StevensOp (line 42)
  - `stevens_operators.py`: function:StevensOp (line 42)
  - `stevens_operators.py`: function:StevensOp (line 44)
  - `stevens_operators.py`: function:StevensOp (line 44)
  - `stevens_operators.py`: function:StevensOp (line 51)
  - `stevens_operators.py`: function:StevensOp (line 51)
  - `stevens_operators.py`: function:StevensOp (line 53)
  - `stevens_operators.py`: function:StevensOp (line 53)
  - `stevens_operators.py`: function:StevensOp (line 55)
  - `stevens_operators.py`: function:StevensOp (line 55)
  - `stevens_operators.py`: function:StevensOp (line 55)
  - `stevens_operators.py`: function:StevensOp (line 55)
  - `stevens_operators.py`: function:StevensOp (line 57)
  - `stevens_operators.py`: function:StevensOp (line 57)
  - `stevens_operators.py`: function:StevensOp (line 61)
  - `stevens_operators.py`: function:StevensOp (line 61)
  - `stevens_operators.py`: function:StevensOp (line 63)
  - `stevens_operators.py`: function:StevensOp (line 63)
  - `stevens_operators.py`: function:StevensOp (line 65)
  - `stevens_operators.py`: function:StevensOp (line 65)
  - `stevens_operators.py`: function:StevensOp (line 65)
  - `stevens_operators.py`: function:StevensOp (line 65)
  - `stevens_operators.py`: function:StevensOp (line 70)
  - `stevens_operators.py`: function:StevensOp (line 70)
  - `stevens_operators.py`: function:StevensOp (line 72)
  - `stevens_operators.py`: function:StevensOp (line 72)
  - `stevens_operators.py`: function:StevensOp (line 74)
  - `stevens_operators.py`: function:StevensOp (line 74)
  - `stevens_operators.py`: function:StevensOp (line 74)
  - `stevens_operators.py`: function:StevensOp (line 74)
  - `stevens_operators.py`: function:StevensOp (line 76)
  - `stevens_operators.py`: function:StevensOp (line 76)
  - `stevens_operators.py`: function:StevensOp (line 77)
  - `stevens_operators.py`: function:StevensOp (line 77)
  - `stevens_operators.py`: function:StevensOp (line 79)
  - `stevens_operators.py`: function:StevensOp (line 79)
  - `stevens_operators.py`: function:StevensOp (line 79)
  - `stevens_operators.py`: function:StevensOp (line 80)
  - `stevens_operators.py`: function:StevensOp (line 80)
  - `stevens_operators.py`: function:StevensOp (line 80)
  - `stevens_operators.py`: function:StevensOp (line 82)
  - `stevens_operators.py`: function:StevensOp (line 82)
  - `stevens_operators.py`: function:StevensOp (line 82)
  - `stevens_operators.py`: function:StevensOp (line 87)
  - `stevens_operators.py`: function:StevensOp (line 87)
  - `stevens_operators.py`: function:StevensOp (line 89)
  - `stevens_operators.py`: function:StevensOp (line 89)
  - `stevens_operators.py`: function:StevensOp (line 91)
  - `stevens_operators.py`: function:StevensOp (line 91)
  - `stevens_operators.py`: function:StevensOp (line 91)
  - `stevens_operators.py`: function:StevensOp (line 91)
  - `stevens_operators.py`: function:StevensOp (line 93)
  - `stevens_operators.py`: function:StevensOp (line 93)
  - `stevens_operators.py`: function:StevensOp (line 94)
  - `stevens_operators.py`: function:StevensOp (line 94)
  - `stevens_operators.py`: function:StevensOp (line 96)
  - `stevens_operators.py`: function:StevensOp (line 96)
  - `stevens_operators.py`: function:StevensOp (line 96)
  - `stevens_operators.py`: function:StevensOp (line 97)
  - `stevens_operators.py`: function:StevensOp (line 97)
  - `stevens_operators.py`: function:StevensOp (line 97)
  </details>

- `Jplus()` *(line 14)* - 5 usages
  <details><summary>View all 5 usages</summary>

  - `operators.py`: class:Ket → method:Jx (line 23)
  - `operators.py`: class:Ket → method:Jy (line 26)
  - `operators.py`: class:Operator → method:Jx (line 120)
  - `operators.py`: class:Operator → method:Jy (line 126)
  - `stevens_operators.py`: function:StevensOp (line 9)
  </details>

- `Jminus()` *(line 18)* - 5 usages
  <details><summary>View all 5 usages</summary>

  - `operators.py`: class:Ket → method:Jx (line 23)
  - `operators.py`: class:Ket → method:Jy (line 26)
  - `operators.py`: class:Operator → method:Jx (line 121)
  - `operators.py`: class:Operator → method:Jy (line 127)
  - `stevens_operators.py`: function:StevensOp (line 10)
  </details>

- `Jx()` *(line 22)* - 30 usages
  <details><summary>View all 30 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:__init__ (line 30)
  - `cf_levels.py`: class:CFLevels → method:Hamiltonian (line 59)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 291)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 291)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 377)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 395)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 406)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 424)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 518)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 518)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 530)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 530)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 565)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 575)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 576)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 577)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 578)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 592)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 592)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 604)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 621)
  - `cf_levels.py`: class:OpticalTransition → method:__init__ (line 730)
  - `cf_levels.py`: class:OpticalTransition → method:transition_strength (line 749)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 822)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1290)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1290)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1313)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1314)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1315)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1316)
  </details>

- `Jy()` *(line 25)* - 30 usages
  <details><summary>View all 30 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:__init__ (line 31)
  - `cf_levels.py`: class:CFLevels → method:Hamiltonian (line 59)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 292)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 292)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 378)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 396)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 406)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 425)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 519)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 519)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 531)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 531)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 566)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 580)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 581)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 582)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 583)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 593)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 593)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 604)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 622)
  - `cf_levels.py`: class:OpticalTransition → method:__init__ (line 731)
  - `cf_levels.py`: class:OpticalTransition → method:transition_strength (line 750)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 823)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1290)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1290)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1318)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1319)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1320)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1321)
  </details>

- `R()` *(line 28)* - 2 usages
  <details><summary>View usages</summary>

  - `thermo_functions.py`: function:Cp1T (line 11)
  - `thermo_functions.py`: function:Cp1T (line 16)
  </details>

- `_Rz()` *(line 31)* - 1 usages
  <details><summary>View usages</summary>

  - `operators.py`: class:Ket → method:R (line 29)
  </details>

- `_Ry()` *(line 37)* - 0 usages
- `_WignersFormula()` *(line 46)* - 1 usages
  <details><summary>View usages</summary>

  - `operators.py`: class:Ket → method:_Ry (line 43)
  </details>

- `__mul__()` *(line 64)* - 0 usages
- `__add__()` *(line 71)* - 0 usages

<details><summary>📍 View all 19 class usages</summary>

- `cf_levels.py`: import-statement (line 10)
- `cf_levels.py`: class:CFLevels → method:gsExpectation (line 376)
- `cf_levels.py`: class:CFLevels → method:gsExpectation (line 376)
- `cf_levels.py`: class:CFLevels → method:magnetization (line 421)
- `cf_levels.py`: class:CFLevels → method:magnetization (line 421)
- `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 516)
- `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 516)
- `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 526)
- `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 526)
- `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 620)
- `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 620)
- `cf_levels.py`: class:CFLevels → method:testEigenvectors (line 691)
- `cf_levels.py`: class:CFLevels → method:testEigenvectors (line 691)
- `cf_levels.py`: class:CFLevels → method:testEigenvectors (line 691)
- `cf_levels.py`: class:CFLevels → method:testEigenvectors (line 691)
- `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum (line 920)
- `cf_levels.py`: class:LS_CFLevels → method:neutronSpectrum (line 920)
- `cf_levels.py`: class:LS_CFLevels → method:normalizedNeutronSpectrum (line 967)
- `cf_levels.py`: class:LS_CFLevels → method:normalizedNeutronSpectrum (line 967)
</details>

#### **Operator** *(line 85)*

- **Used by:** 14 files


**Methods:**

- `__init__()` *(line 86)* - 0 usages
- `Jz()` *(line 92)* - 105 usages
  <details><summary>View all 105 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:__init__ (line 32)
  - `cf_levels.py`: class:CFLevels → method:Hamiltonian (line 59)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 293)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 293)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 379)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 397)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 406)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 426)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 520)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 520)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 532)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 532)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 567)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 569)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 570)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 571)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 572)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 594)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 594)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 604)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 623)
  - `cf_levels.py`: class:OpticalTransition → method:__init__ (line 732)
  - `cf_levels.py`: class:OpticalTransition → method:transition_strength (line 751)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 824)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1290)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1290)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1308)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1309)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1310)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1311)
  - `stevens_operators.py`: function:StevensOp (line 8)
  - `stevens_operators.py`: function:StevensOp (line 8)
  - `stevens_operators.py`: function:StevensOp (line 16)
  - `stevens_operators.py`: function:StevensOp (line 25)
  - `stevens_operators.py`: function:StevensOp (line 25)
  - `stevens_operators.py`: function:StevensOp (line 27)
  - `stevens_operators.py`: function:StevensOp (line 29)
  - `stevens_operators.py`: function:StevensOp (line 29)
  - `stevens_operators.py`: function:StevensOp (line 36)
  - `stevens_operators.py`: function:StevensOp (line 36)
  - `stevens_operators.py`: function:StevensOp (line 38)
  - `stevens_operators.py`: function:StevensOp (line 38)
  - `stevens_operators.py`: function:StevensOp (line 40)
  - `stevens_operators.py`: function:StevensOp (line 40)
  - `stevens_operators.py`: function:StevensOp (line 42)
  - `stevens_operators.py`: function:StevensOp (line 42)
  - `stevens_operators.py`: function:StevensOp (line 44)
  - `stevens_operators.py`: function:StevensOp (line 44)
  - `stevens_operators.py`: function:StevensOp (line 51)
  - `stevens_operators.py`: function:StevensOp (line 51)
  - `stevens_operators.py`: function:StevensOp (line 53)
  - `stevens_operators.py`: function:StevensOp (line 53)
  - `stevens_operators.py`: function:StevensOp (line 55)
  - `stevens_operators.py`: function:StevensOp (line 55)
  - `stevens_operators.py`: function:StevensOp (line 55)
  - `stevens_operators.py`: function:StevensOp (line 55)
  - `stevens_operators.py`: function:StevensOp (line 57)
  - `stevens_operators.py`: function:StevensOp (line 57)
  - `stevens_operators.py`: function:StevensOp (line 61)
  - `stevens_operators.py`: function:StevensOp (line 61)
  - `stevens_operators.py`: function:StevensOp (line 63)
  - `stevens_operators.py`: function:StevensOp (line 63)
  - `stevens_operators.py`: function:StevensOp (line 65)
  - `stevens_operators.py`: function:StevensOp (line 65)
  - `stevens_operators.py`: function:StevensOp (line 65)
  - `stevens_operators.py`: function:StevensOp (line 65)
  - `stevens_operators.py`: function:StevensOp (line 70)
  - `stevens_operators.py`: function:StevensOp (line 70)
  - `stevens_operators.py`: function:StevensOp (line 72)
  - `stevens_operators.py`: function:StevensOp (line 72)
  - `stevens_operators.py`: function:StevensOp (line 74)
  - `stevens_operators.py`: function:StevensOp (line 74)
  - `stevens_operators.py`: function:StevensOp (line 74)
  - `stevens_operators.py`: function:StevensOp (line 74)
  - `stevens_operators.py`: function:StevensOp (line 76)
  - `stevens_operators.py`: function:StevensOp (line 76)
  - `stevens_operators.py`: function:StevensOp (line 77)
  - `stevens_operators.py`: function:StevensOp (line 77)
  - `stevens_operators.py`: function:StevensOp (line 79)
  - `stevens_operators.py`: function:StevensOp (line 79)
  - `stevens_operators.py`: function:StevensOp (line 79)
  - `stevens_operators.py`: function:StevensOp (line 80)
  - `stevens_operators.py`: function:StevensOp (line 80)
  - `stevens_operators.py`: function:StevensOp (line 80)
  - `stevens_operators.py`: function:StevensOp (line 82)
  - `stevens_operators.py`: function:StevensOp (line 82)
  - `stevens_operators.py`: function:StevensOp (line 82)
  - `stevens_operators.py`: function:StevensOp (line 87)
  - `stevens_operators.py`: function:StevensOp (line 87)
  - `stevens_operators.py`: function:StevensOp (line 89)
  - `stevens_operators.py`: function:StevensOp (line 89)
  - `stevens_operators.py`: function:StevensOp (line 91)
  - `stevens_operators.py`: function:StevensOp (line 91)
  - `stevens_operators.py`: function:StevensOp (line 91)
  - `stevens_operators.py`: function:StevensOp (line 91)
  - `stevens_operators.py`: function:StevensOp (line 93)
  - `stevens_operators.py`: function:StevensOp (line 93)
  - `stevens_operators.py`: function:StevensOp (line 94)
  - `stevens_operators.py`: function:StevensOp (line 94)
  - `stevens_operators.py`: function:StevensOp (line 96)
  - `stevens_operators.py`: function:StevensOp (line 96)
  - `stevens_operators.py`: function:StevensOp (line 96)
  - `stevens_operators.py`: function:StevensOp (line 97)
  - `stevens_operators.py`: function:StevensOp (line 97)
  - `stevens_operators.py`: function:StevensOp (line 97)
  </details>

- `Jplus()` *(line 101)* - 5 usages
  <details><summary>View all 5 usages</summary>

  - `operators.py`: class:Ket → method:Jx (line 23)
  - `operators.py`: class:Ket → method:Jy (line 26)
  - `operators.py`: class:Operator → method:Jx (line 120)
  - `operators.py`: class:Operator → method:Jy (line 126)
  - `stevens_operators.py`: function:StevensOp (line 9)
  </details>

- `Jminus()` *(line 110)* - 5 usages
  <details><summary>View all 5 usages</summary>

  - `operators.py`: class:Ket → method:Jx (line 23)
  - `operators.py`: class:Ket → method:Jy (line 26)
  - `operators.py`: class:Operator → method:Jx (line 121)
  - `operators.py`: class:Operator → method:Jy (line 127)
  - `stevens_operators.py`: function:StevensOp (line 10)
  </details>

- `Jx()` *(line 119)* - 30 usages
  <details><summary>View all 30 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:__init__ (line 30)
  - `cf_levels.py`: class:CFLevels → method:Hamiltonian (line 59)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 291)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 291)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 377)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 395)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 406)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 424)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 518)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 518)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 530)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 530)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 565)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 575)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 576)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 577)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 578)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 592)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 592)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 604)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 621)
  - `cf_levels.py`: class:OpticalTransition → method:__init__ (line 730)
  - `cf_levels.py`: class:OpticalTransition → method:transition_strength (line 749)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 822)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1290)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1290)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1313)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1314)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1315)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1316)
  </details>

- `Jy()` *(line 125)* - 30 usages
  <details><summary>View all 30 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:__init__ (line 31)
  - `cf_levels.py`: class:CFLevels → method:Hamiltonian (line 59)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 292)
  - `cf_levels.py`: class:CFLevels → method:_transition (line 292)
  - `cf_levels.py`: class:CFLevels → method:gsExpectation (line 378)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 396)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 406)
  - `cf_levels.py`: class:CFLevels → method:magnetization (line 425)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 519)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 519)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 531)
  - `cf_levels.py`: class:CFLevels → method:susceptibilityPert (line 531)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 566)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 580)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 581)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 582)
  - `cf_levels.py`: class:CFLevels → method:gtensor (line 583)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 593)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 593)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 604)
  - `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 622)
  - `cf_levels.py`: class:OpticalTransition → method:__init__ (line 731)
  - `cf_levels.py`: class:OpticalTransition → method:transition_strength (line 750)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 823)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1290)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1290)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1318)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1319)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1320)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensor (line 1321)
  </details>

- `__add__()` *(line 130)* - 0 usages
- `__radd__()` *(line 138)* - 0 usages
- `__sub__()` *(line 146)* - 0 usages
- `__mul__()` *(line 154)* - 0 usages
- `__rmul__()` *(line 162)* - 0 usages
- `__pow__()` *(line 170)* - 0 usages
- `__neg__()` *(line 177)* - 0 usages
- `__repr__()` *(line 182)* - 0 usages

<details><summary>📍 View all 14 class usages</summary>

- `cf_levels.py`: import-statement (line 10)
- `cf_levels.py`: class:CFLevels → method:__init__ (line 30)
- `cf_levels.py`: class:CFLevels → method:__init__ (line 31)
- `cf_levels.py`: class:CFLevels → method:__init__ (line 32)
- `cf_levels.py`: class:CFLevels → method:Hamiltonian (line 59)
- `cf_levels.py`: class:CFLevels → method:Hamiltonian (line 59)
- `cf_levels.py`: class:CFLevels → method:Hamiltonian (line 59)
- `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 592)
- `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 593)
- `cf_levels.py`: class:CFLevels → method:gtensorzeeman (line 594)
- `stevens_operators.py`: import-statement (line 2)
- `stevens_operators.py`: function:StevensOp (line 8)
- `stevens_operators.py`: function:StevensOp (line 9)
- `stevens_operators.py`: function:StevensOp (line 10)
</details>

#### **LSOperator** *(line 205)*

- **Used by:** 25 files


**Methods:**

- `__init__()` *(line 207)* - 0 usages
- `Lz()` *(line 217)* - 17 usages
  <details><summary>View all 17 usages</summary>

  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 811)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 811)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 813)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 824)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 827)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1352)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1352)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1356)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1357)
  - `constants.py`: function:PFgamma (line 271)
  - `constants.py`: function:PFgamma (line 275)
  - `constants.py`: function:PFgamma (line 276)
  - `constants.py`: function:PFgamma (line 277)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 233)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 233)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 235)
  - `operators.py`: class:LSOperator → method:magnetization (line 363)
  </details>

- `Lplus()` *(line 226)* - 2 usages
  <details><summary>View usages</summary>

  - `operators.py`: class:LSOperator → method:Lx (line 245)
  - `operators.py`: class:LSOperator → method:Ly (line 251)
  </details>

- `Lminus()` *(line 235)* - 2 usages
  <details><summary>View usages</summary>

  - `operators.py`: class:LSOperator → method:Lx (line 246)
  - `operators.py`: class:LSOperator → method:Ly (line 252)
  </details>

- `Lx()` *(line 244)* - 13 usages
  <details><summary>View all 13 usages</summary>

  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 809)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 809)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 813)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 822)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 825)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1350)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1350)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1356)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1357)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 231)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 231)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 235)
  - `operators.py`: class:LSOperator → method:magnetization (line 361)
  </details>

- `Ly()` *(line 250)* - 13 usages
  <details><summary>View all 13 usages</summary>

  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 810)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 810)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 813)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 823)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 826)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1351)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1351)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1356)
  - `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1357)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 232)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 232)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 235)
  - `operators.py`: class:LSOperator → method:magnetization (line 362)
  </details>

- `Sz()` *(line 258)* - 9 usages
  <details><summary>View all 9 usages</summary>

  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 808)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 808)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 813)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 824)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 827)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 230)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 230)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 235)
  - `operators.py`: class:LSOperator → method:magnetization (line 366)
  </details>

- `Splus()` *(line 267)* - 2 usages
  <details><summary>View usages</summary>

  - `operators.py`: class:LSOperator → method:Sx (line 286)
  - `operators.py`: class:LSOperator → method:Sy (line 292)
  </details>

- `Sminus()` *(line 276)* - 2 usages
  <details><summary>View usages</summary>

  - `operators.py`: class:LSOperator → method:Sx (line 287)
  - `operators.py`: class:LSOperator → method:Sy (line 293)
  </details>

- `Sx()` *(line 285)* - 9 usages
  <details><summary>View all 9 usages</summary>

  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 806)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 806)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 813)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 822)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 825)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 228)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 228)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 235)
  - `operators.py`: class:LSOperator → method:magnetization (line 364)
  </details>

- `Sy()` *(line 291)* - 9 usages
  <details><summary>View all 9 usages</summary>

  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 807)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 807)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 813)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 823)
  - `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 826)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 229)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 229)
  - `ligands.py`: class:LS_Ligands → method:__init__ (line 235)
  - `operators.py`: class:LSOperator → method:magnetization (line 365)
  </details>

- `__add__()` *(line 297)* - 0 usages
- `__radd__()` *(line 305)* - 0 usages
- `__sub__()` *(line 313)* - 0 usages
- `__mul__()` *(line 321)* - 0 usages
- `__rmul__()` *(line 329)* - 0 usages
- `__pow__()` *(line 337)* - 0 usages
- `__neg__()` *(line 344)* - 0 usages
- `__repr__()` *(line 349)* - 0 usages
- `magnetization()` *(line 355)* - 48 usages
  <details><summary>View all 48 usages</summary>

  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 459)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 460)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 461)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 462)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 468)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 469)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 470)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 471)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 477)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 478)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 479)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 480)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 488)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 489)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 490)
  - `cf_levels.py`: class:CFLevels → method:susceptibility (line 491)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1122)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1123)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1124)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1125)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1131)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1132)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1133)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1134)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1140)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1141)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1142)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1143)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1151)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1152)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1153)
  - `cf_levels.py`: class:LS_CFLevels → method:susceptibility (line 1154)
  - `operators.py`: class:LSOperator → method:susceptibility (line 423)
  - `operators.py`: class:LSOperator → method:susceptibility (line 424)
  - `operators.py`: class:LSOperator → method:susceptibility (line 425)
  - `operators.py`: class:LSOperator → method:susceptibility (line 426)
  - `operators.py`: class:LSOperator → method:susceptibility (line 432)
  - `operators.py`: class:LSOperator → method:susceptibility (line 433)
  - `operators.py`: class:LSOperator → method:susceptibility (line 434)
  - `operators.py`: class:LSOperator → method:susceptibility (line 435)
  - `operators.py`: class:LSOperator → method:susceptibility (line 441)
  - `operators.py`: class:LSOperator → method:susceptibility (line 442)
  - `operators.py`: class:LSOperator → method:susceptibility (line 443)
  - `operators.py`: class:LSOperator → method:susceptibility (line 444)
  - `operators.py`: class:LSOperator → method:susceptibility (line 452)
  - `operators.py`: class:LSOperator → method:susceptibility (line 453)
  - `operators.py`: class:LSOperator → method:susceptibility (line 454)
  - `operators.py`: class:LSOperator → method:susceptibility (line 455)
  </details>

- `susceptibility()` *(line 413)* - 0 usages

<details><summary>📍 View all 25 class usages</summary>

- `cf_levels.py`: import-statement (line 10)
- `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 795)
- `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 795)
- `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 806)
- `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 807)
- `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 808)
- `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 809)
- `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 810)
- `cf_levels.py`: class:LS_CFLevels → method:__init__ (line 811)
- `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1350)
- `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1351)
- `cf_levels.py`: class:LS_CFLevels → method:gtensorperturb (line 1352)
- `ligands.py`: import-statement (line 10)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 228)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 229)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 230)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 231)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 232)
- `ligands.py`: class:LS_Ligands → method:__init__ (line 233)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 339)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 339)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 415)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 415)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 491)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 491)
</details>

---

## 📁 plot_ligands.py {#plot_ligands-py}

### 🏛️ Classes

#### **atomplot** *(line 53)*

- **Used by:** 0 file


**Methods:**

- `__init__()` *(line 54)* - 0 usages
- `plotatoms()` *(line 63)* - 1 usages
  <details><summary>View usages</summary>

  - `plot_ligands.py`: function:plotPCF (line 10)
  </details>

- `plotaxes()` *(line 87)* - 1 usages
  <details><summary>View usages</summary>

  - `plot_ligands.py`: function:plotPCF (line 11)
  </details>

- `plotabc()` *(line 99)* - 0 usages
- `_flatten()` *(line 114)* - 12 usages
  <details><summary>View all 12 usages</summary>

  - `plot_ligands.py`: class:atomplot → method:plotaxes (line 90)
  - `plot_ligands.py`: class:atomplot → method:plotaxes (line 91)
  - `plot_ligands.py`: class:atomplot → method:plotaxes (line 92)
  - `plot_ligands.py`: class:atomplot → method:plotaxes (line 95)
  - `plot_ligands.py`: class:atomplot → method:plotaxes (line 96)
  - `plot_ligands.py`: class:atomplot → method:plotaxes (line 97)
  - `plot_ligands.py`: class:atomplot → method:plotabc (line 105)
  - `plot_ligands.py`: class:atomplot → method:plotabc (line 106)
  - `plot_ligands.py`: class:atomplot → method:plotabc (line 107)
  - `plot_ligands.py`: class:atomplot → method:plotabc (line 110)
  - `plot_ligands.py`: class:atomplot → method:plotabc (line 111)
  - `plot_ligands.py`: class:atomplot → method:plotabc (line 112)
  </details>

### 📊 Functions

#### **plotPCF()** *(line 5)*

- **Used by:** 3 files

<details><summary>View all 3 usages</summary>

- `cif_symmetry_import.py`: import-statement (line 3)
- `cif_symmetry_import.py`: function:FindPointGroupSymOps (line 319)
- `cif_symmetry_import.py`: function:FindPointGroupSymOps (line 319)
</details>

#### **exportLigandCif()** *(line 15)*

- **Used by:** 3 files

<details><summary>View all 3 usages</summary>

- `ligands.py`: import-statement (line 5)
- `ligands.py`: class:Ligands → method:exportCif (line 68)
- `ligands.py`: class:Ligands → method:exportCif (line 68)
</details>

---

## 📁 rescale_CEF.py {#rescale_CEF-py}

### 📊 Functions

#### **rescaleCEF()** *(line 3)*

- **Used by:** 0 file

---

## 📁 stevens_operators.py {#stevens_operators-py}

### 📊 Functions

#### **StevensOp()** *(line 5)*

- **Used by:** 14 files

<details><summary>View all 14 usages</summary>

- `cf_levels.py`: import-statement (line 11)
- `cf_levels.py`: class:CFLevels → method:Bdict (line 47)
- `cf_levels.py`: class:CFLevels → method:Bdict (line 47)
- `ligands.py`: import-statement (line 7)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 131)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 131)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 136)
- `ligands.py`: class:Ligands → method:PointChargeModel (line 136)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 333)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 333)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 409)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 409)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 485)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 485)
</details>

#### **LS_StevensOp()** *(line 102)*

- **Used by:** 10 files

<details><summary>View all 10 usages</summary>

- `cf_levels.py`: import-statement (line 11)
- `cf_levels.py`: class:LS_CFLevels → method:Bdict (line 840)
- `cf_levels.py`: class:LS_CFLevels → method:Bdict (line 840)
- `ligands.py`: import-statement (line 7)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 330)
- `ligands.py`: class:LS_Ligands → method:PointChargeModel (line 330)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 406)
- `ligands.py`: class:LS_Ligands → method:TMPointChargeModel (line 406)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 482)
- `ligands.py`: class:LS_Ligands → method:UnknownTMPointChargeModel (line 482)
</details>

---

## 📁 thermo_functions.py {#thermo_functions-py}

### 📊 Functions

#### **partition_func()** *(line 5)*

- **Used by:** 0 file

#### **Cp_from_CEF()** *(line 9)*

- **Used by:** 0 file

#### **Cp1T()** *(line 10)*

- **Used by:** 0 file

---

## 📁 wybourne_stevens.py {#wybourne_stevens-py}

### 📊 Functions

#### **WybourneToStevens()** *(line 3)*

- **Used by:** 0 file

#### **StevensToWybourne()** *(line 14)*

- **Used by:** 0 file

---

## 📁 undefined_files\PCF_misc_functions.py {#undefined_files-PCF_misc_functions-py}

### 📊 Functions

#### **backgroundfunction()** *(line 7)*

- **Used by:** 0 file

#### **importfile()** *(line 32)*

- **Used by:** 0 file

#### **importGridfile()** *(line 39)*

- **Used by:** 0 file

---

## 📈 Summary Statistics

- **Total files analyzed:** 20
- **Total entities discovered:** 208
- **Total cross-file dependencies:** 1067

### 🏆 Most Frequently Used Entities

| Entity | Type | File | Usage Count |
|--------|------|------|-------------|
| `Jz` | METHOD | operators.py | 105 |
| `Jz` | METHOD | operators.py | 105 |
| `magnetization` | METHOD | cf_levels.py | 48 |
| `magnetization` | METHOD | cf_levels.py | 48 |
| `magnetization` | METHOD | operators.py | 48 |
| `Jx` | METHOD | operators.py | 30 |
| `Jx` | METHOD | operators.py | 30 |
| `Jy` | METHOD | operators.py | 30 |
| `Jy` | METHOD | operators.py | 30 |
| `k_B` | VARIABLE | constants.py | 29 |

### 🔗 Files with Most Dependencies

| File | Dependency Count |
|------|------------------|
| cf_levels.py | 516 |
| ligands.py | 164 |
| stevens_operators.py | 158 |
| operators.py | 88 |
| import_CIF.py | 34 |

### ⚠️ Potentially Unused Entities (34 total)

<details><summary>View all unused entities</summary>

| Entity | Type | File | Line |
|--------|------|------|------|
| `_spec` | VARIABLE | cf_levels.py | 703 |
| `OpticalTransition` | CLASS | cf_levels.py | 710 |
| `LandeGFactor` | FUNCTION | cf_levels.py | 765 |
| `findRotationAxis` | FUNCTION | cif_symmetry_import.py | 332 |
| `makeSymOpMatrix` | FUNCTION | cif_symmetry_import.py | 376 |
| `ION_NUMS_TRANS_METAL` | CONSTANT | constants.py | 19 |
| `_tesseral_dispatch` | FUNCTION | constants.py | 117 |
| `PFgamma` | FUNCTION | constants.py | 240 |
| `SPIN_ORBIT_COUPLING_CM` | CONSTANT | constants.py | 341 |
| `resultfunc` | FUNCTION | create_fit_function.py | 51 |
| `makeCurveFitFunction` | FUNCTION | create_fit_function.py | 61 |
| `resultfunc` | FUNCTION | create_fit_function.py | 106 |
| `importRE_FF` | FUNCTION | form_factors.py | 4 |
| `importCIF` | FUNCTION | import_CIF.py | 12 |
| `checkTMexist` | FUNCTION | import_CIF.py | 118 |
| `printLaTexCEFparams` | FUNCTION | latex_cef_print.py | 3 |
| `MomIntertia` | FUNCTION | moments_of_inertia.py | 5 |
| `selectZaxisMI` | FUNCTION | moments_of_inertia.py | 21 |
| `ContinuousShapeMeasure` | FUNCTION | moments_of_inertia.py | 40 |
| `anglesToVector` | FUNCTION | moments_of_inertia.py | 52 |
| `rotationMatrix` | FUNCTION | moments_of_inertia.py | 56 |
| `rotateArbAxis` | FUNCTION | moments_of_inertia.py | 72 |
| `findZaxis_SOM_rotation` | FUNCTION | moments_of_inertia.py | 80 |
| `fitfun` | FUNCTION | moments_of_inertia.py | 83 |
| `atomplot` | CLASS | plot_ligands.py | 53 |
| `rescaleCEF` | FUNCTION | rescale_CEF.py | 3 |
| `partition_func` | FUNCTION | thermo_functions.py | 5 |
| `Cp_from_CEF` | FUNCTION | thermo_functions.py | 9 |
| `Cp1T` | FUNCTION | thermo_functions.py | 10 |
| `WybourneToStevens` | FUNCTION | wybourne_stevens.py | 3 |
| `StevensToWybourne` | FUNCTION | wybourne_stevens.py | 14 |
| `backgroundfunction` | FUNCTION | undefined_files\PCF_misc_functions.py | 7 |
| `importfile` | FUNCTION | undefined_files\PCF_misc_functions.py | 32 |
| `importGridfile` | FUNCTION | undefined_files\PCF_misc_functions.py | 39 |

</details>