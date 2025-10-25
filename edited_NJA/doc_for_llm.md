# NJA-CFS Technical Reference

## Overview
Crystal field software for electronic structure calculations of d^n and f^n configurations. Handles interelectronic repulsion, spin-orbit coupling, ligand/crystal field effects, and magnetic properties.

## Core Classes

### `Wigner_coeff`
Angular momentum coupling coefficients and rotation matrices.

**Methods:**
- `threej_symbol(matrix)` → float: 3j-symbol via Racah formula
- `sixj_symbol(matrix)` → float: 6j-symbol via Racah formula  
- `Wigner_Dmatrix(l, m1, m, alpha, beta, gamma)` → complex: Rotation matrix element (Euler angles, radians)
- `Wigner_Dmatrix_quat_complete(l, R, bin, dict, coeff)` → ndarray: Full Wigner D-matrix from quaternion

**Input:** Angular momentum quantum numbers, rotation parameters
**Output:** Coupling coefficients, rotation matrices

### `CFP`
Coefficients of fractional parentage for electron configurations.

**Init:** `CFP(conf, dic_cfp, dic_LS_inv_almost=None)`
- `conf`: String like 'd3', 'f5'
- Handles almost-closed shells automatically

**Methods:**
- `cfp(v, L, S, name)` → ndarray: CFP values for state (v, L, S)

### `RME`
Reduced matrix elements for tensor operators.

**Init:** `RME(state, conf, dic_cfp, labels, dic_LS, dic_LS_inv_almost)`
- `state`: [v, L, S, v1, L1, S1]

**Methods:**
- `Uk(k, ...)` → float: Reduced matrix element for rank-k operator
- `V1k(k=1)` → float: Spin-dependent matrix element

### `Hamiltonian`
Matrix elements for different interaction terms.

**Init:** `Hamiltonian(state, labels, conf, dic_cfp, tables, dic_LS, dic_LS_almost)`
- `state`: [v, L, S, v1, L1, S1, J, M, J1, M1]

**Methods:**
- `electrostatic_int(basis, F0, F2, F4, F6, evaluation, tab_ee)` → float: Electron repulsion (cm⁻¹)
- `SO_coupling(zeta, k, evaluation)` → float: Spin-orbit coupling (cm⁻¹)
- `LF_contribution(dic_kq, evaluation)` → complex: Crystal/ligand field term (cm⁻¹)
- `Zeeman(field, k, evaluation, MM)` → complex: Zeeman splitting (cm⁻¹)

**Crystal Field Input:** `dic_bkq = {'2': {'0': value, '1': value, ...}, '4': {...}, ...}`

### `calculation`
Main computation engine - builds and diagonalizes Hamiltonian.

**Init:** `calculation(conf, ground_only=False, TAB=False, wordy=True)`
- `conf`: 'd1'-'d9', 'f1'-'f13'
- `ground_only`: Restricts to ground multiplet
- `TAB`: Use tabulated RMEs (faster) vs explicit calculation

**Attributes:**
- `basis`: Complete basis set [2S, L, 2J, 2M, seniority]
- `dic_LS`: State labels dictionary
- `basis_l`, `basis_l_JM`: Human-readable labels

**Methods:**
- `MatrixH(elem, F0, F2, F4, F6, zeta, k, dic_V, dic_bkq, dic_AOM, PCM, field, cfp_angles, ...)` → ndarray
  - `elem`: List from ['Hee', 'Hso', 'Hcf', 'Hz']
  - Returns: Eigenvalues (row 0) and eigenvectors (subsequent rows)
  
- `build_matrix(elem, ...)` → ndarray: Constructs Hamiltonian matrix
  
- `reduce_basis(conf, roots, dic, contributes)`: Reduces basis by spin multiplicity
  
- `ground_state_calc(ground)`: Extracts ground multiplet basis

**Output:** Energy levels (cm⁻¹), wavefunctions, projections

### `Magnetics`
Magnetic properties from calculated states.

**Init:** `Magnetics(calc, contributes, par, wordy=False)`
- `calc`: calculation object
- `contributes`: ['Hee', 'Hso', 'Hcf']
- `par`: Parameter dictionary

**Methods:**
- `mag_moment(k, evaluation)` → ndarray: Magnetic moment matrix [3, N, N]
- `calc_LS(k)` → tuple: L and S operator matrices
- `effGval(levels, v_matrix)` → tuple: Effective g-values and axes for Kramers doublet
- `susceptibility_field(fields, temp, delta)` → tuple: χ(T) and M(T) along field directions
- `susceptibility_B_ord1(fields, temp, basis, LF_matrix, delta)` → tuple: χ-tensor via Richardson extrapolation

**Inputs:**
- `fields`: Nx3 array (Tesla)
- `temp`: Temperature (K)
- `delta`: Differentiation step (T)

**Outputs:** SI units (m³, A·m²)

## Key Functions

### Crystal Field Conversions

**`from_Vint_to_Bkq(dic_V, conf)`**
- Converts one-electron integrals to Wybourne Bkq parameters
- Input: `dic_V = {'11': value, '21': value, '22': value, ...}` (cm⁻¹)
- Output: `dic_bkq = {'2': {'0': ..., '1': ..., '-1': ...}, ...}` (cm⁻¹)

**`from_AOM_to_Vint(dic_AOM, conf)`**
- Angular Overlap Model → one-electron integrals
- Input: `dic_AOM = {'ligand_name': [e_σ, e_π_s, e_π_c, θ, φ, χ]}` (angles in degrees)

**`calc_Bkq(data, conf, sph_flag, sth_param)`**
- Point charge model → Bkq
- Input: `data` = Nx5 array [label, x, y, z, charge] or spherical coords
- `sth_param`: Apply Sternheimer shielding

**`from_Aqkrk_to_Bkq(Aqkrk, revers)`**
- Stevens ↔ Wybourne conversion

**`rota_LF(l, dic_Bkq, A, B, C)`**
- Rotate Bkq by Euler angles (radians, ZYZ convention)

**`rota_LF_quat(l, dic_Bkq, R, dict, coeff)`**
- Rotate Bkq by quaternion [w, x, y, z]

### Basis and Projections

**`Full_basis(conf)`** → tuple
- Complete basis set for configuration
- Returns: (basis array, dic_LS, basis_l, basis_l_JM)

**`projection_basis(basis2, labels, bin, J_label)`** → dict
- Projects eigenstates onto free-ion basis
- Returns: `{state_index: {label: percentage, ...}}`

**`terms_labels(conf)`** → list
- Free-ion term symbols (Nielson-Koster order)

**`terms_basis(conf)`** → list
- [2S, L, seniority] for each term

### Free Ion Parameters

**`free_ion_param_f(conf)`** → dict
- Slater-Condon F^k and ζ for f^2-f^12 (cm⁻¹)
- Source: Goerller-Walrand & Binnemans

**`free_ion_param_f_HF(conf)`** → dict
- Hartree-Fock parameters f^1-f^13
- Source: Ma et al. 2014

**`free_ion_param_AB(conf)`** → dict
- Parameters from Abragam & Bleaney

### File I/O

**`read_AILFT_orca6(filename, conf, method, return_V, rotangle_V, return_orcamatrix)`**
- Extracts AILFT matrices from ORCA 6 output
- `method`: 'CASSCF', 'NEVPT2', etc.
- `rotangle_V`: Optional rotation [α, β, γ] or quaternion
- Returns: Parameter dictionary or V-matrix

**`cfp_from_file(conf)`** → dict
- Reads CFP from tables/cfp_d_conf.txt or cfp_f_conf.txt

**`read_matrix_from_file(conf, closed_shell)`** → dict
- Loads tabulated RMEs (Uk, V1k)

**`read_data(filename, sph_flag)`** → ndarray
- Ligand coordinates and charges from file

### Coordinate Transformations

**`from_car_to_sph(coord)`** → ndarray
- Cartesian → spherical [r, θ, φ] (radians)

**`from_sph_to_car(coord_sph)`** → ndarray
- Spherical → Cartesian

**`rotate_vectors(vectors, angle, axis)`** → ndarray
- Rotates Nx3 array around 'x', 'y', or 'z'

**`points_on_sphere(num_pts, figure, angles)`**
- Golden spiral point distribution on unit sphere

### Magnetic Calculations (Numba-accelerated)

**`mag_moment(basis)`** → ndarray
- Constructs μ-matrices [3, N, N] from basis
- Returns: -k·L - ge·S operators

**`M_vector(field_vec, mu_matrix, LF_matrix, basis, temp)`** → ndarray
- Magnetization vector at temperature (Bohr magnetons)

**`susceptibility_B_ord1(fields, temp, basis, LF_matrix, delta)`** → tuple
- χ-tensor and error via Ridders' differentiation

**`add_Zeeman(field_vec, basis, LF_matrix)`** → ndarray
- Adds Zeeman term to LF matrix

**`effGval(levels, mu_matrix, v_matrix)`** → tuple
- g-values for Kramers doublet

**`calc_torque(B0, T, LF_matrix, basis, plane, step, figure, show_fig)`** → tuple
- Magnetic torque τ vs rotation angle

### Utilities

**`state_legend(L_str, inv)`**
- Convert 'S', 'P', 'D', ... ↔ 0, 1, 2, ...

**`almost_closed_shells(name)`** → int
- Maps d^6-d^9, f^8-f^13 to effective electron count

**`ground_term_legend(conf)`** → str
- Ground state term symbol

**`r_expect(k, conf)`** → float
- ⟨r^k⟩ expectation values (atomic units)

**`sigma_k(k, conf)`** → float
- Sternheimer shielding parameters

**`Stev_coeff(k, conf)`** → float
- Stevens coefficients (α, β, γ)

**`conv_Aqkrk_bkq(l, m)`** → float
- Stevens → Wybourne conversion factor

### Plotting

**`plot_energy_levels(eigenvalues, ax, color, label, tolerance, offset, delta, ylabel)`**
- Energy level diagram with degeneracy grouping

**`fig_tensor_rep_1(tensor, n_points)`**
- 3D surface plot of tensor magnitude M(n) = n^T·T·n

**`fig_susc_field(conf, dic_Bkq, temp, n_points, delta)`**
- χ(field direction) visualization

**`level_fig_tot(E_matrix, theories, proj_LS_dict, proj_prev_dict)`**
- Multi-theory correlation diagram

## Typical Workflows

### 1. Ground State Calculation
```python
calc = calculation('f7', ground_only=True, TAB=True)
dic = {'dic_bkq': {...}, 'F2': ..., 'F4': ..., 'F6': ..., 'zeta': ...}
result = calc.MatrixH(['Hee', 'Hso', 'Hcf'], **dic)
```

### 2. Magnetic Properties
```python
Magn = Magnetics(calc, ['Hee', 'Hso', 'Hcf'], dic)
chi_tensor, err = Magn.susceptibility_B_ord1(fields, temp, calc.basis, LF_matrix)
g_values, g_axes = Magn.effGval([0, 1])  # Ground Kramers doublet
```

### 3. From AILFT Data
```python
dic = read_AILFT_orca6('output.out', 'f7', method='CASSCF')
# dic contains F2, F4, F6, zeta, dic_bkq ready for calculation
```

### 4. Point Charge Model
```python
data = read_data('ligands.txt', sph_flag=False)
dic_bkq = calc_Bkq(data, 'f7', sph_flag=False, sth_param=True)
```

## Units
- **Energy:** cm⁻¹
- **Magnetic field:** Tesla
- **Susceptibility:** m³ (SI) or cm³·K/mol (CGS with conversion)
- **Magnetization:** Bohr magnetons
- **Angles:** Radians (unless specified)
- **Distances:** Ångström (file input)

## Configuration Strings
- d-orbitals: 'd1' through 'd9'
- f-orbitals: 'f1' through 'f13'
- Automatically handles almost-closed shells (d^6-d^9, f^8-f^13)