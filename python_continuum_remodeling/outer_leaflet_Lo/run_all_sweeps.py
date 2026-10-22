"""
Master Runner Script for HermesLoouter (Protein on outer leaflet of Lo, bare Ld).
Executes:
1. Smoke test verification.
2. Single-domain expanded sweeps (126 simulations).
3. Two-domain expanded sweeps (126 simulations).
4. High-resolution publication figure generation.
"""

import os
import sys
import time

current_dir = os.path.dirname(os.path.abspath(__file__))
single_dir = os.path.join(current_dir, "singledomain")
two_dir = os.path.join(current_dir, "twodomain")

sys.path.insert(0, current_dir)
sys.path.insert(0, single_dir)
sys.path.insert(0, two_dir)

from smoke_test import run_smoke_test
from param_sweeps import run_all_single_domain_sweeps
from two_domains_sweeps import run_all_two_domain_sweeps

def main():
    t_start = time.time()
    print("=" * 80)
    print("HERMES LO-OUTER PRODUCTION SIMULATION PIPELINE")
    print("=" * 80)
    
    # 1. Run Smoke Test
    print("\n>>> STEP 0: RUNNING PRE-PRODUCTION SMOKE TEST AUDIT...")
    run_smoke_test()
    
    # 2. Single Domain Production Sweep
    print("\n>>> STEP 1: RUNNING SINGLE-DOMAIN SWEEPS (126 RUNS)...")
    run_all_single_domain_sweeps(out_dir=os.path.join(current_dir, "single"))
    
    # 3. Two Domain Production Sweep
    print("\n>>> STEP 2: RUNNING TWO-DOMAIN SWEEPS (126 RUNS)...")
    run_all_two_domain_sweeps(out_dir=os.path.join(current_dir, "double"))
    
    # 4. Generate Figures
    print("\n>>> STEP 3: REGENERATING PUBLICATION FIGURES...")
    try:
        from singledomain.plot_param_sweeps import plot_param_sweeps
        plot_param_sweeps(results_dir=os.path.join(current_dir, "single"))
    except Exception as e:
        print(f"Warning: Single domain plot failed: {e}")
        
    try:
        from twodomain.plot_two_domains import plot_two_domains
        plot_two_domains(results_dir=os.path.join(current_dir, "double"))
    except Exception as e:
        print(f"Warning: Two domain plot failed: {e}")
        
    try:
        from plot_comparison_single_two import plot_comparison
        plot_comparison(root_dir=current_dir)
    except Exception as e:
        print(f"Warning: Comparison plot failed: {e}")
        
    t_tot = time.time() - t_start
    print("\n" + "=" * 80)
    print(f"ALL 252 PRODUCTION SIMULATIONS AND PLOTS COMPLETED IN {t_tot:.1f} SECONDS!")
    print("=" * 80)

if __name__ == "__main__":
    main()
