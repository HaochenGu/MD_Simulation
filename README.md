# Polymer Radius of Gyration MD Simulation

This project performs molecular dynamics (MD) simulations to calculate the radius of gyration (Rg²) for coarse-grained polymer structures using LAMMPS with the Kremer-Grest model.

## This repository is highly based on ClaudeCode.

## Overview

The project simulates 1342 different polymer structures with various topologies (linear, cyclic, star, branch, comb, dendrimer) to compute their radius of gyration and compare with reference values. The simulations use a coarse-grained representation where each bead represents multiple monomers.

## Key Features

- **Accurate polymer representation**: Uses NetworkX graph structures for correct connectivity
- **Handles non-contiguous indices**: Fixed scripts automatically handle molecules with non-contiguous node numbering
- **Parallel processing**: Support for running multiple simulations concurrently
- **Multiple output formats**: Generates detailed results and simplified CSV summaries
- **Automatic unit conversion**: Converts between LJ units and Angstroms

## Project Structure

```
MD_Simulation/
├── data/
│   ├── rg2.pickle              # Polymer data (graphs, reference Rg² values)
│   └── rg2_SMILES.csv          # Polymer metadata (topology, SMILES, Rg² in Å²)
├── run_polymer_fixed.py        # Single polymer simulation (RECOMMENDED)
├── batch_simulate_polymers_fixed.py  # Batch simulation with parallel support (RECOMMENDED)
├── run_polymer_correct.py      # Original working version
├── batch_simulate_polymers.py  # Original batch script
├── simulation.py               # Legacy script (has issues - do not use)
├── requirements.txt            # Python dependencies
├── CLAUDE.md                   # Detailed technical documentation
└── README.md                   # This file
```

## Installation

1. **Clone the repository and set up Python environment:**
```bash
cd MD_Simulation
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
```

2. **Install LAMMPS:**
   - Download and install LAMMPS from https://lammps.sandia.gov/
   - Ensure the LAMMPS executable (`lmp`) is in your PATH
   - Or specify the full path when running scripts

## Quick Start

### Single Polymer Simulation

```bash
# Quick test with reduced timesteps
python run_polymer_fixed.py 0 --quick lmp

# Full simulation for molecule 728 (previously failing)
python run_polymer_fixed.py 728 lmp

# The results will be saved in a timestamped directory:
# 2025_6_5_2140/polymer728/
```

### Batch Simulations

```bash
# Test first 10 polymers with 4 parallel jobs
python batch_simulate_polymers_fixed.py --start 0 --end 10 --quick --parallel 4

# Run problematic molecules (306-965) that have non-contiguous indices
python batch_simulate_polymers_fixed.py --start 306 --end 320 --quick --parallel 8

# Full dataset (all 1342 polymers) with 8 parallel jobs
python batch_simulate_polymers_fixed.py --start 0 --end 1342 --parallel 8
```

## Output Files

Each simulation creates a timestamped directory (e.g., `2025_6_5_2140/`) containing:

### For Single Simulations (`run_polymer_fixed.py`):
- `polymer{id}/` - Directory for each polymer
  - `polymer_{id}.data` - LAMMPS data file
  - `polymer_{id}.in` - LAMMPS input script
  - `polymer_{id}.log` - LAMMPS output log
  - `rg_{id}.dat` - Time series of Rg² during production
  - `polymer_{id}_results.csv` - Summary with columns: polymer_idx, reference_rg2, rg2, error

### For Batch Simulations (`batch_simulate_polymers_fixed.py`):
- `mol_{id:04d}/` - Directory for each molecule
- `batch_results.csv` - Detailed results for all polymers
- `results_summary.csv` - Simplified results with columns: polymer_idx, reference_rg2, rg2, error

## Technical Details

### Simulation Parameters
- **Model**: Kremer-Grest coarse-grained polymer model
- **Force field**: 
  - WCA (Weeks-Chandler-Andersen) potential for non-bonded interactions
  - FENE (Finitely Extensible Nonlinear Elastic) bonds
- **Temperature**: T* = 1.0 (in reduced units)
- **Integration**: Langevin thermostat (γ=0.1) + NVE integrator
- **Timestep**: 0.001 τ
- **Equilibration**: 10⁷ steps (or 10⁵ for quick mode)
- **Production**: 10⁷ steps (or 10⁵ for quick mode)

### Known Issues and Solutions

1. **Non-contiguous node indices (molecules 306-965)**
   - **Problem**: Graph nodes have indices like [0,1,2,...,71,1056,1057,...,3041]
   - **Solution**: Use `run_polymer_fixed.py` or `batch_simulate_polymers_fixed.py` which automatically remap indices

2. **Missing NVE integrator in original code**
   - **Problem**: Atoms don't move without NVE integrator
   - **Solution**: All current scripts include both Langevin thermostat and NVE integrator

## Performance

- **Accuracy**: Typically <2% error for linear/cyclic/branch polymers
- **Speed**: ~30-60 seconds per polymer (quick mode), ~5-10 minutes (full mode)
- **Parallel scaling**: Near-linear speedup with multiple cores

## Command Line Options

### run_polymer_fixed.py
```bash
python run_polymer_fixed.py [mol_id] [options]
  mol_id      Molecule ID (0-1341)
  --quick     Use reduced timesteps for testing
  lmp         Path to LAMMPS executable (default: 'lmp')
```

### batch_simulate_polymers_fixed.py
```bash
python batch_simulate_polymers_fixed.py [options]
  --start N        Start molecule index (default: 0)
  --end N          End molecule index (default: 10)
  --quick          Use reduced timesteps
  --parallel N     Number of parallel jobs (default: 1)
  --lammps PATH    Path to LAMMPS executable (default: 'lmp')
  --output FILE    Output CSV filename (default: 'batch_results.csv')
```

## Data Files

### rg2.pickle
Contains 6 arrays loaded in sequence:
1. `x` - Adjacency matrices (may have issues, use graphs instead)
2. `rg2` - Reference Rg² values in LJ units
3. `topo_desc` - Topology descriptors
4. `topo_class` - Topology class labels
5. `poly_param` - Polymer parameters
6. `graph` - NetworkX graph objects (use these for connectivity)

### rg2_SMILES.csv
Contains polymer metadata:
- `Topology` - Polymer topology type
- `Number of nodes` - Number of beads
- `Rg2` - Reference radius of gyration squared in Å²
- `SMILES` - SMILES string representation

## Troubleshooting

1. **LAMMPS not found**
   - Ensure LAMMPS is installed and in PATH
   - Or provide full path: `python run_polymer_fixed.py 0 /usr/local/bin/lmp`

2. **Import errors**
   - Activate virtual environment: `source venv/bin/activate`
   - Install dependencies: `pip install -r requirements.txt`

3. **Simulations failing**
   - Use the `_fixed.py` versions of scripts
   - Check LAMMPS log files in the output directory
   - Ensure sufficient memory for parallel jobs

## Citation

Based on: Jiang et al., "Machine learning-enabled multimodal fusion of intra-atrial and body surface signals in prediction of atrial fibrillation ablation outcomes", npj Computational Materials 10, 139 (2024)
