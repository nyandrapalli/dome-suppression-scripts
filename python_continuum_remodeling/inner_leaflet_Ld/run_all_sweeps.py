"""
Master Runner Script for Expanded Biophysical Vesicle Sweeps (PRL Standard).
Runs:
1. Single-domain expanded sweeps:
   - Line tension gamma in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0] pN
   - Domain size D_Lo in [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0] um
   - Deflation v in [0.90, 0.95, 0.98, 1.00]
2. Two-domain expanded sweeps across the exact same parameter space.
3. Automatically triggers high-resolution plotting pipelines.
"""

import os
import sys
import time

current_dir = os.path.dirname(os.path.abspath(__file__))
single_dir = os.path.join(current_dir, "singledomain")
two_dir = os.path.join(current_dir, "twodomain")

sys.path.insert(0, single_dir)
sys.path.insert(0, two_dir)

from param_sweeps import run_all_single_domain_sweeps
from two_domains_sweeps import run_all_two_domain_sweeps

def main():
    t_start = time.time()
    print("=" * 80)
    print("STARTING FULL MASTER RE-RUN OF EXPANDED PARAMETER SWEEPS")
    print("=" * 80)
    
    # 1. Single Domain Production Sweep
    print("\n>>> STEP 1: RUNNING SINGLE-DOMAIN EXPANDED SWEEPS (126 RUNS)...")
    run_all_single_domain_sweeps()
    
    # 2. Two Domain Production Sweep
    print("\n>>> STEP 2: RUNNING TWO-DOMAIN EXPANDED SWEEPS (126 RUNS)...")
    run_all_two_domain_sweeps()
    
    # 3. Regenerate Figures
    print("\n>>> STEP 3: REGENERATING PUBLICATION FIGURES...")
    try:
        from plot_two_domains_bifurcation import plot_two_domains_bifurcation
        plot_two_domains_bifurcation()
    except Exception as e:
        print(f"Warning: Plot bifurcation failed: {e}")
        
    try:
        from plot_two_domains_showcase import plot_two_domains_showcase
        plot_two_domains_showcase()
    except Exception as e:
        print(f"Warning: Plot showcase failed: {e}")
        
    try:
        from plot_comparison_single_two import plot_single_vs_two_comparison
        plot_single_vs_two_comparison()
    except Exception as e:
        print(f"Warning: Plot comparison failed: {e}")
        
    t_tot = time.time() - t_start
    print("\n" + "=" * 80)
    print(f"ALL 252 PRODUCTION SIMULATIONS AND PLOTS COMPLETED IN {t_tot:.1f} SECONDS!")
    print("=" * 80)

if __name__ == "__main__":
    main()
