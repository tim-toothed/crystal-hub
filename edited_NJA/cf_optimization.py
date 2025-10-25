"""
Crystal Field Level Optimization for HoFeO₃
==========================================

This script performs systematic optimization of crystal field parameters
to match experimental energy levels from literature.

Process: Crystal Field Parameter Fitting (Phenomenological Analysis)
- Adjusts model parameters while preserving structural symmetry
- Standard procedure in rare-earth spectroscopy
- Accounts for limitations of point charge model
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution
import nja_cfs_red as nja
from cif_parser import parse_cif_to_nja
from typing import Dict, List, Tuple, Optional
import json
from datetime import datetime

# Conversion factor
CM_TO_MEV = 0.123984193

# ============================================================================
# EXPERIMENTAL DATA (Target values from Ovsyanikov 2020)
# ============================================================================

EXPERIMENTAL_LEVELS_MEV = np.array([
    0.81,   # Ground state (taking their ~0.81 meV as reference)
    10.62,
    15.64,
    21.67,
    27.65,
    34.65,
    40.59,
    44.62
])

# We can also include the full 17 levels if available
# For now using the clearly visible ones from the figure

class CFOptimizer:
    """
    Crystal Field Parameter Optimizer
    
    Performs systematic fitting of CF parameters to experimental data.
    This is a standard procedure in rare-earth spectroscopy called
    "phenomenological crystal field analysis".
    """
    
    def __init__(self, cif_file: str, metal_element: str = 'Ho'):
        self.cif_file = cif_file
        self.metal_element = metal_element
        self.configuration = 'f10'  # Ho3+
        
        # Store base calculation parameters
        self.base_coordination_data = None
        self.base_free_ion_params = None
        
        # Results storage
        self.optimization_history = []
        self.best_params = None
        self.best_score = np.inf
        
    def calculate_levels(self, 
                        cutoff: float,
                        charge_scaling: float = 1.0,
                        beta: float = 1.0,
                        f2_scaling: float = 1.0,
                        zeta_scaling: float = 1.0) -> np.ndarray:
        """
        Calculate CF levels with given parameters.
        
        Parameters:
        -----------
        cutoff : float
            Coordination cutoff distance (Å)
        charge_scaling : float
            Scaling factor for effective charges (typically 0.8-1.2)
        beta : float
            Nephelauxetic ratio (accounts for covalency, typically 0.90-1.00)
        f2_scaling : float
            Scaling for F² parameter (typically 0.90-1.00)
        zeta_scaling : float
            Scaling for spin-orbit coupling (typically 0.95-1.05)
            
        Returns:
        --------
        energy_levels : ndarray
            Calculated energy levels in meV
        """
        
        # Parse CIF with specified cutoff
        coordination_data = parse_cif_to_nja(
            cif_file=self.cif_file,
            metal_center=self.metal_element,
            coordination_cutoff=cutoff,
            charge_dict=None,
            output_file=None,
            auto_detect=False
        )
        
        # Scale charges if requested
        if charge_scaling != 1.0:
            scaled_coords = []
            for coord in coordination_data:
                scaled_coord = list(coord)
                scaled_coord[4] = coord[4] * charge_scaling  # Scale the charge
                scaled_coords.append(tuple(scaled_coord))
            coordination_data = scaled_coords

        # Convert to numpy array before calculating B_kq
        coordination_data = np.array(coordination_data)

        # Calculate B_kq parameters
        dic_bkq = nja.calc_Bkq(
            coordination_data,
            self.configuration,
            sph_flag=False,
            sth_param=True
        )
        
        # Get free ion parameters
        free_ion_params = nja.free_ion_param_f_HF(self.configuration)
        
        # Apply scaling factors
        free_ion_params['F2'] *= f2_scaling * beta
        free_ion_params['F4'] *= beta
        free_ion_params['F6'] *= beta
        free_ion_params['zeta'] *= zeta_scaling
        
        # Calculate energy levels
        calc = nja.calculation(
            self.configuration,
            ground_only=True,
            TAB=True,
            wordy=False
        )
        
        result = calc.MatrixH(
            ['Hee', 'Hso', 'Hcf'],
            F2=free_ion_params['F2'],
            F4=free_ion_params['F4'],
            F6=free_ion_params['F6'],
            zeta=free_ion_params['zeta'],
            dic_bkq=dic_bkq,
            wordy=False
        )
        
        if isinstance(result, tuple):
            eigenvalues = result[0][0, :].real
        else:
            eigenvalues = result[0, :].real
        
        # Convert to meV relative to ground state
        energies_relative_cm = eigenvalues - eigenvalues[0]
        energies_relative_meV = energies_relative_cm * CM_TO_MEV
        
        return energies_relative_meV
    
    def cost_function(self, params: np.ndarray) -> float:
        """
        Cost function for optimization.
        
        Measures how well calculated levels match experimental data.
        Uses Hungarian algorithm concept - matches each calculated level
        to closest experimental level.
        
        Parameters:
        -----------
        params : ndarray
            [cutoff, charge_scaling, beta, f2_scaling, zeta_scaling]
            
        Returns:
        --------
        cost : float
            Root mean square deviation (RMSD) between calculated and experimental levels
        """
        
        cutoff, charge_scaling, beta, f2_scaling, zeta_scaling = params
        
        try:
            # Calculate levels
            calc_levels = self.calculate_levels(
                cutoff, charge_scaling, beta, f2_scaling, zeta_scaling
            )
            
            # Match calculated to experimental levels
            # For each experimental level, find closest calculated level
            min_distances = []
            used_calc_indices = set()
            
            for exp_level in EXPERIMENTAL_LEVELS_MEV:
                min_dist = np.inf
                best_idx = -1
                
                for i, calc_level in enumerate(calc_levels):
                    if i in used_calc_indices:
                        continue
                    dist = abs(calc_level - exp_level)
                    if dist < min_dist:
                        min_dist = dist
                        best_idx = i
                
                if best_idx != -1:
                    min_distances.append(min_dist)
                    used_calc_indices.add(best_idx)
                else:
                    # Penalty if we can't match all experimental levels
                    min_distances.append(10.0)  # 10 meV penalty
            
            # RMSD
            rmsd = np.sqrt(np.mean(np.array(min_distances)**2))
            
            # Penalty for having too many extra levels
            n_extra = len(calc_levels) - len(EXPERIMENTAL_LEVELS_MEV)
            if n_extra > 0:
                rmsd += 0.5 * n_extra  # Penalty for extra levels
            
            # Store in history
            self.optimization_history.append({
                'params': params.copy(),
                'rmsd': rmsd,
                'calc_levels': calc_levels.copy()
            })
            
            if rmsd < self.best_score:
                self.best_score = rmsd
                self.best_params = params.copy()
                print(f"New best RMSD: {rmsd:.4f} meV")
                print(f"  Cutoff: {cutoff:.3f} Å")
                print(f"  Charge scaling: {charge_scaling:.3f}")
                print(f"  β: {beta:.4f}")
                print(f"  F² scaling: {f2_scaling:.4f}")
                print(f"  ζ scaling: {zeta_scaling:.4f}")
            
            return rmsd
            
        except Exception as e:
            print(f"Error in calculation: {e}")
            return 1000.0  # Large penalty for failed calculations
    
    def optimize_grid_search(self, 
                            cutoff_range: Tuple[float, float] = (2.5, 3.5),
                            n_points: int = 20) -> Dict:
        """
        Systematic grid search over cutoff distance.
        
        This is often the most important parameter to adjust.
        
        Parameters:
        -----------
        cutoff_range : tuple
            (min_cutoff, max_cutoff) in Angstroms
        n_points : int
            Number of points to sample
            
        Returns:
        --------
        results : dict
            Dictionary with cutoff values and corresponding RMSDs
        """
        
        print("=" * 70)
        print("GRID SEARCH: Cutoff Distance Optimization")
        print("=" * 70)
        
        cutoffs = np.linspace(cutoff_range[0], cutoff_range[1], n_points)
        rmsd_values = []
        
        for cutoff in cutoffs:
            # Use default values for other parameters
            params = np.array([cutoff, 1.0, 1.0, 1.0, 1.0])
            rmsd = self.cost_function(params)
            rmsd_values.append(rmsd)
            
        results = {
            'cutoffs': cutoffs,
            'rmsds': np.array(rmsd_values),
            'best_cutoff': cutoffs[np.argmin(rmsd_values)],
            'best_rmsd': np.min(rmsd_values)
        }
        
        # Plot results
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(cutoffs, rmsd_values, 'bo-', linewidth=2, markersize=8)
        ax.axvline(results['best_cutoff'], color='r', linestyle='--', 
                   label=f"Best cutoff: {results['best_cutoff']:.3f} Å")
        ax.set_xlabel('Coordination Cutoff (Å)', fontsize=12)
        ax.set_ylabel('RMSD (meV)', fontsize=12)
        ax.set_title('Effect of Coordination Cutoff on CF Level Agreement', fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.legend()
        plt.tight_layout()
        plt.savefig('cutoff_optimization.png', dpi=300)
        plt.show()
        
        return results
    
    def optimize_full(self, 
                     initial_cutoff: float = 2.89,
                     method: str = 'differential_evolution') -> Dict:
        """
        Full multi-parameter optimization.
        
        Uses either:
        - 'differential_evolution': Global optimizer, slower but more robust
        - 'nelder-mead': Local optimizer, faster but may get stuck
        
        Parameters:
        -----------
        initial_cutoff : float
            Starting cutoff distance
        method : str
            Optimization method
            
        Returns:
        --------
        results : dict
            Optimized parameters and final RMSD
        """
        
        print("=" * 70)
        print(f"FULL OPTIMIZATION: {method}")
        print("=" * 70)
        
        # Parameter bounds: [cutoff, charge_scaling, beta, f2_scaling, zeta_scaling]
        bounds = [
            (2.3, 3.5),      # cutoff: reasonable range for Ho-O distances
            (0.7, 1.3),      # charge_scaling: effective charges can vary
            (0.90, 1.00),    # beta: nephelauxetic ratio (covalency reduction)
            (0.90, 1.00),    # F2_scaling: typical reduction due to environment
            (0.95, 1.05)     # zeta_scaling: small adjustments to spin-orbit
        ]
        
        if method == 'differential_evolution':
            # Global optimization - slower but thorough
            result = differential_evolution(
                self.cost_function,
                bounds,
                maxiter=100,
                popsize=15,
                tol=0.01,
                seed=42,
                workers=1,
                disp=True
            )
            optimal_params = result.x
            final_rmsd = result.fun
            
        else:  # nelder-mead or other scipy methods
            # Local optimization - faster
            x0 = np.array([initial_cutoff, 1.0, 0.98, 0.95, 1.0])
            result = minimize(
                self.cost_function,
                x0,
                method='Nelder-Mead',
                options={'maxiter': 1000, 'xatol': 0.001, 'fatol': 0.01}
            )
            optimal_params = result.x
            final_rmsd = result.fun
        
        # Calculate final levels with optimal parameters
        final_levels = self.calculate_levels(*optimal_params)
        
        results = {
            'optimal_cutoff': optimal_params[0],
            'optimal_charge_scaling': optimal_params[1],
            'optimal_beta': optimal_params[2],
            'optimal_f2_scaling': optimal_params[3],
            'optimal_zeta_scaling': optimal_params[4],
            'final_rmsd': final_rmsd,
            'calculated_levels': final_levels,
            'experimental_levels': EXPERIMENTAL_LEVELS_MEV
        }
        
        return results
    
    def plot_comparison(self, results: Dict, save_path: str = 'level_comparison.png'):
        """
        Plot comparison between calculated and experimental levels.
        """
        
        calc_levels = results['calculated_levels']
        exp_levels = results['experimental_levels']
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 8))
        
        # Left panel: Energy level diagram
        ax1.hlines(exp_levels, 0.3, 0.5, colors='red', linewidth=3, 
                   label='Experimental (Ovsyanikov 2020)')
        ax1.hlines(calc_levels, 0.6, 0.8, colors='blue', linewidth=3, 
                   label='Calculated (Optimized)')
        
        ax1.set_xlim(0, 1)
        ax1.set_ylabel('Energy (meV)', fontsize=12)
        ax1.set_title('CF Levels: Experimental vs Calculated', fontsize=14)
        ax1.legend(loc='upper right')
        ax1.set_xticks([])
        ax1.grid(axis='y', alpha=0.3)
        
        # Right panel: Correlation plot
        # Match closest levels
        matched_pairs = []
        used_calc = set()
        for exp in exp_levels:
            diffs = [abs(exp - calc) for i, calc in enumerate(calc_levels) if i not in used_calc]
            if diffs:
                min_idx = np.argmin([abs(exp - calc) for calc in calc_levels])
                matched_pairs.append((exp, calc_levels[min_idx]))
                used_calc.add(min_idx)
        
        if matched_pairs:
            exp_matched, calc_matched = zip(*matched_pairs)
            ax2.scatter(exp_matched, calc_matched, s=100, alpha=0.6)
            
            # Ideal line
            max_val = max(max(exp_matched), max(calc_matched))
            ax2.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, label='Ideal agreement')
            
            ax2.set_xlabel('Experimental Energy (meV)', fontsize=12)
            ax2.set_ylabel('Calculated Energy (meV)', fontsize=12)
            ax2.set_title(f"Level Correlation (RMSD = {results['final_rmsd']:.3f} meV)", fontsize=14)
            ax2.legend()
            ax2.grid(True, alpha=0.3)
            ax2.set_aspect('equal')
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"\nComparison plot saved to: {save_path}")
    
    def save_results(self, results: Dict, filename: str = 'optimization_results.json'):
        """
        Save optimization results to file.
        """
        
        # Convert numpy arrays to lists for JSON serialization
        serializable_results = {
            'timestamp': datetime.now().isoformat(),
            'compound': f'{self.metal_element}FeO3',
            'configuration': self.configuration,
            'optimal_parameters': {
                'cutoff_angstrom': float(results['optimal_cutoff']),
                'charge_scaling': float(results['optimal_charge_scaling']),
                'beta_nephelauxetic': float(results['optimal_beta']),
                'F2_scaling': float(results['optimal_f2_scaling']),
                'zeta_scaling': float(results['optimal_zeta_scaling'])
            },
            'final_rmsd_meV': float(results['final_rmsd']),
            'calculated_levels_meV': results['calculated_levels'].tolist(),
            'experimental_levels_meV': results['experimental_levels'].tolist(),
            'level_comparison': []
        }
        
        # Add level-by-level comparison
        for i, (exp, calc) in enumerate(zip(
            results['experimental_levels'], 
            results['calculated_levels'][:len(results['experimental_levels'])]
        )):
            serializable_results['level_comparison'].append({
                'level': i,
                'experimental_meV': float(exp),
                'calculated_meV': float(calc),
                'difference_meV': float(calc - exp)
            })
        
        with open(filename, 'w') as f:
            json.dump(serializable_results, f, indent=2)
        
        print(f"\nResults saved to: {filename}")
    
    def print_summary(self, results: Dict):
        """
        Print detailed summary of optimization results.
        """
        
        print("\n" + "=" * 70)
        print("OPTIMIZATION SUMMARY")
        print("=" * 70)
        
        print(f"\nCompound: {self.metal_element}FeO₃")
        print(f"Configuration: {self.configuration} ({self.metal_element}³⁺)")
        
        print("\n" + "-" * 70)
        print("OPTIMAL PARAMETERS")
        print("-" * 70)
        print(f"  Coordination cutoff:     {results['optimal_cutoff']:.4f} Å")
        print(f"  Effective charge scaling: {results['optimal_charge_scaling']:.4f}")
        print(f"  β (nephelauxetic ratio): {results['optimal_beta']:.4f}")
        print(f"  F² scaling factor:       {results['optimal_f2_scaling']:.4f}")
        print(f"  ζ scaling factor:        {results['optimal_zeta_scaling']:.4f}")
        
        print("\n" + "-" * 70)
        print("AGREEMENT WITH EXPERIMENT")
        print("-" * 70)
        print(f"  Final RMSD: {results['final_rmsd']:.4f} meV")
        
        print("\n" + "-" * 70)
        print("LEVEL-BY-LEVEL COMPARISON (meV)")
        print("-" * 70)
        print(f"{'Level':<8} {'Experimental':<15} {'Calculated':<15} {'Δ (Calc-Exp)':<15}")
        print("-" * 70)
        
        calc_levels = results['calculated_levels']
        exp_levels = results['experimental_levels']
        
        # Match and compare
        used_calc = set()
        for i, exp in enumerate(exp_levels):
            # Find closest calculated level
            min_diff = np.inf
            best_calc = None
            best_idx = -1
            
            for j, calc in enumerate(calc_levels):
                if j not in used_calc:
                    diff = abs(calc - exp)
                    if diff < min_diff:
                        min_diff = diff
                        best_calc = calc
                        best_idx = j
            
            if best_calc is not None:
                used_calc.add(best_idx)
                diff = best_calc - exp
                print(f"{i+1:<8} {exp:<15.3f} {best_calc:<15.3f} {diff:<+15.3f}")
            else:
                print(f"{i+1:<8} {exp:<15.3f} {'---':<15} {'---':<15}")
        
        # Print any unmatched calculated levels
        unmatched = [calc for i, calc in enumerate(calc_levels) if i not in used_calc]
        if unmatched:
            print("\nUnmatched calculated levels:")
            for calc in unmatched:
                print(f"         {'---':<15} {calc:<15.3f}")
        
        print("\n" + "=" * 70)


# ============================================================================
# MAIN OPTIMIZATION WORKFLOW
# ============================================================================

def main():
    """
    Main optimization workflow for HoFeO₃ crystal field levels.
    """
    
    # Configuration
    CIF_FILE = r'C:\Users\Timur\Documents\Python_Scripts\crystal-hub\cif_files\0.994_HoFeO3.mcif'
    
    print("=" * 70)
    print("CRYSTAL FIELD PARAMETER OPTIMIZATION FOR HoFeO₃")
    print("=" * 70)
    print("\nThis script performs 'phenomenological crystal field analysis'")
    print("to optimize agreement between calculated and experimental CF levels.")
    print("\nProcess:")
    print("1. Grid search over coordination cutoff")
    print("2. Multi-parameter optimization")
    print("3. Comparison with experimental data")
    print("=" * 70)
    
    # Initialize optimizer
    optimizer = CFOptimizer(CIF_FILE, metal_element='Ho')
    
    # Step 1: Grid search over cutoff (usually most important)
    print("\n" + "=" * 70)
    print("STEP 1: Grid Search - Finding optimal coordination cutoff")
    print("=" * 70)
    
    grid_results = optimizer.optimize_grid_search(
        cutoff_range=(2.5, 3.2),
        n_points=30
    )
    
    print(f"\nBest cutoff from grid search: {grid_results['best_cutoff']:.3f} Å")
    print(f"RMSD at best cutoff: {grid_results['best_rmsd']:.3f} meV")
    
    # Step 2: Full multi-parameter optimization
    print("\n" + "=" * 70)
    print("STEP 2: Full Multi-Parameter Optimization")
    print("=" * 70)
    
    # Use best cutoff from grid search as starting point
    optimization_results = optimizer.optimize_full(
        initial_cutoff=grid_results['best_cutoff'],
        method='differential_evolution'  # Try 'nelder-mead' for faster (but less thorough) optimization
    )
    
    # Step 3: Analyze and visualize results
    print("\n" + "=" * 70)
    print("STEP 3: Analysis and Visualization")
    print("=" * 70)
    
    optimizer.print_summary(optimization_results)
    optimizer.plot_comparison(optimization_results)
    optimizer.save_results(optimization_results)
    
    print("\n" + "=" * 70)
    print("OPTIMIZATION COMPLETE!")
    print("=" * 70)
    print("\nNext steps:")
    print("1. Check 'level_comparison.png' for visual agreement")
    print("2. Review 'optimization_results.json' for detailed parameters")
    print("3. If RMSD is still high, you can:")
    print("   - Adjust parameter bounds in optimize_full()")
    print("   - Run longer optimization (increase maxiter)")
    print("   - Try different optimization methods")
    print("   - Consider that perfect agreement may not be possible")
    print("     (the point charge model has fundamental limitations)")
    

if __name__ == '__main__':
    main()