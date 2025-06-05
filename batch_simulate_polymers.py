#!/usr/bin/env python3
"""
Batch simulation script for all polymers using the correct approach.
Uses graph structures from pickle file for accurate polymer representation.
"""

import numpy as np
import pickle
import pandas as pd
import networkx as nx
import os
import subprocess
import time
from datetime import datetime
import argparse
from multiprocessing import Pool
import sys

# Add the directory to path so we can import the simulation function
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from run_polymer_correct import run_polymer_simulation


def process_polymer(args):
    """Wrapper function for multiprocessing."""
    mol_id, lammps_exe, quick, work_dir, data_dir = args
    
    # Change to work directory
    original_dir = os.getcwd()
    mol_dir = os.path.join(work_dir, f'mol_{mol_id:04d}')
    os.makedirs(mol_dir, exist_ok=True)
    os.chdir(mol_dir)
    
    try:
        result = run_polymer_simulation(mol_id, lammps_exe, quick, data_dir)
    except Exception as e:
        print(f"Error processing molecule {mol_id}: {e}")
        result = {'mol_id': mol_id, 'success': False, 'error': str(e)}
    finally:
        os.chdir(original_dir)
    
    return result


def main():
    parser = argparse.ArgumentParser(description='Batch MD simulations for polymers')
    parser.add_argument('--start', type=int, default=0, help='Start molecule index')
    parser.add_argument('--end', type=int, default=10, help='End molecule index')
    parser.add_argument('--lammps', default='lmp', help='LAMMPS executable')
    parser.add_argument('--quick', action='store_true', help='Use shorter simulations')
    parser.add_argument('--parallel', type=int, default=1, help='Number of parallel jobs')
    parser.add_argument('--output', default='batch_results.csv', help='Output CSV file')
    
    args = parser.parse_args()
    
    # Create work directory with format YYYY_M_D_HHMM
    now = datetime.now()
    timestamp = f"{now.year}_{now.month}_{now.day}_{now.hour:02d}{now.minute:02d}"
    work_dir = timestamp
    os.makedirs(work_dir, exist_ok=True)
    
    print(f"Starting batch simulation in {work_dir}")
    print(f"Processing molecules {args.start} to {args.end-1}")
    print(f"Quick mode: {args.quick}")
    print(f"Parallel jobs: {args.parallel}")
    
    # Get absolute path to data directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_dir, 'data')
    
    # Prepare job list
    jobs = []
    for mol_id in range(args.start, args.end):
        jobs.append((mol_id, args.lammps, args.quick, work_dir, data_dir))
    
    # Run simulations
    results = []
    start_time = time.time()
    
    if args.parallel > 1:
        # Parallel execution
        with Pool(processes=args.parallel) as pool:
            for i, result in enumerate(pool.imap(process_polymer, jobs)):
                results.append(result)
                if result.get('success'):
                    print(f"[{i+1}/{len(jobs)}] Mol {result['mol_id']}: "
                          f"{result['topology']}, Error = {result['error_percent']:.1f}%")
                else:
                    print(f"[{i+1}/{len(jobs)}] Mol {result['mol_id']}: FAILED")
                
                # Save intermediate results
                if (i + 1) % 10 == 0:
                    df = pd.DataFrame(results)
                    df.to_csv(os.path.join(work_dir, args.output), index=False)
    else:
        # Serial execution
        for i, job in enumerate(jobs):
            result = process_polymer(job)
            results.append(result)
            
            if result.get('success'):
                print(f"[{i+1}/{len(jobs)}] Mol {result['mol_id']}: "
                      f"{result['topology']}, Error = {result['error_percent']:.1f}%")
            else:
                print(f"[{i+1}/{len(jobs)}] Mol {result['mol_id']}: FAILED")
            
            # Save intermediate results
            if (i + 1) % 10 == 0:
                df = pd.DataFrame(results)
                df.to_csv(os.path.join(work_dir, args.output), index=False)
    
    # Save final results
    df_results = pd.DataFrame(results)
    df_results.to_csv(os.path.join(work_dir, args.output), index=False)
    
    # Create simplified CSV with requested columns
    simplified_results = []
    for r in results:
        if r.get('success'):
            simplified_results.append({
                'polymer_idx': r['mol_id'],
                'reference_rg2': r['rg2_ref_ang'],
                'rg2': r['rg2_sim'] * r['sigma']**2,  # Convert to Angstroms²
                'error': r['error_percent']
            })
    
    df_simplified = pd.DataFrame(simplified_results)
    simplified_csv = os.path.join(work_dir, 'results_summary.csv')
    df_simplified.to_csv(simplified_csv, index=False)
    print(f"\nSimplified results saved to: {simplified_csv}")
    
    # Summary statistics
    total_time = time.time() - start_time
    successful = df_results[df_results['success'] == True]
    
    print(f"\n{'='*60}")
    print("BATCH SIMULATION SUMMARY")
    print(f"{'='*60}")
    print(f"Total molecules: {len(df_results)}")
    print(f"Successful: {len(successful)}")
    print(f"Failed: {len(df_results) - len(successful)}")
    print(f"Total time: {total_time/60:.1f} minutes")
    print(f"Average time per molecule: {total_time/len(df_results):.1f} seconds")
    
    if len(successful) > 0:
        print(f"\nError statistics:")
        print(f"  Mean error: {successful['error_percent'].mean():.2f}%")
        print(f"  Median error: {successful['error_percent'].median():.2f}%")
        print(f"  Min error: {successful['error_percent'].min():.2f}%")
        print(f"  Max error: {successful['error_percent'].max():.2f}%")
        
        # Group by topology
        print(f"\nError by topology:")
        for topo in sorted(successful['topology'].unique()):
            topo_data = successful[successful['topology'] == topo]
            print(f"  {topo:10s}: {topo_data['error_percent'].mean():6.2f}% "
                  f"(n={len(topo_data):3d}, min={topo_data['error_percent'].min():.1f}%, "
                  f"max={topo_data['error_percent'].max():.1f}%)")
        
        # Show conversion factors
        print(f"\nConversion factors (σ in Å) by topology:")
        for topo in sorted(successful['topology'].unique()):
            topo_data = successful[successful['topology'] == topo]
            print(f"  {topo:10s}: {topo_data['sigma'].mean():.3f} ± {topo_data['sigma'].std():.3f} Å")
    
    print(f"\nResults saved to: {os.path.join(work_dir, args.output)}")


if __name__ == '__main__':
    main()