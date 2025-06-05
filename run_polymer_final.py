#!/usr/bin/env python3
"""
Final robust implementation for polymer simulation.
Handles all edge cases including disconnected graphs and perfect rings.
"""

import numpy as np
import pickle
import pandas as pd
import networkx as nx
import os
import subprocess
import time
import sys
from datetime import datetime


def fix_graph_indices(g):
    """Remap graph nodes to use contiguous indices 0 to n-1."""
    old_nodes = sorted(g.nodes())
    node_map = {old: new for new, old in enumerate(old_nodes)}
    
    # Create new graph with remapped nodes
    g_new = nx.Graph()
    for i in range(len(old_nodes)):
        g_new.add_node(i)
    
    # Add edges with remapped indices
    for u, v in g.edges():
        g_new.add_edge(node_map[u], node_map[v])
    
    return g_new


def place_ring(n_atoms, target_bond_length=0.9):
    """Place atoms in a ring with specified bond length."""
    # For a regular n-gon with side length L, radius R = L/(2*sin(pi/n))
    radius = target_bond_length / (2 * np.sin(np.pi / n_atoms))
    coords = np.zeros((n_atoms, 3))
    
    for i in range(n_atoms):
        angle = 2 * np.pi * i / n_atoms
        coords[i] = [radius * np.cos(angle), radius * np.sin(angle), 0.0]
    
    return coords


def place_component_bfs(graph, component, target_bond_length=0.9):
    """Place a connected component using BFS."""
    coords = {}
    placed = set()
    
    # Start from node with highest degree
    degrees = {n: graph.degree(n) for n in component}
    start = max(degrees, key=degrees.get)
    
    coords[start] = np.array([0.0, 0.0, 0.0])
    placed.add(start)
    
    queue = [start]
    while queue:
        current = queue.pop(0)
        current_pos = coords[current]
        
        # Get unplaced neighbors
        neighbors = [n for n in graph.neighbors(current) if n in component and n not in placed]
        
        for i, neighbor in enumerate(neighbors):
            # Calculate position
            if len(neighbors) == 1:
                theta = np.random.uniform(0, 2*np.pi)
                phi = np.random.uniform(np.pi/3, 2*np.pi/3)
            else:
                theta = 2 * np.pi * i / len(neighbors)
                phi = np.pi/2
            
            # Add small perturbation
            theta += np.random.uniform(-0.1, 0.1)
            phi += np.random.uniform(-0.1, 0.1)
            
            dx = target_bond_length * np.sin(phi) * np.cos(theta)
            dy = target_bond_length * np.sin(phi) * np.sin(theta)
            dz = target_bond_length * np.cos(phi)
            
            coords[neighbor] = current_pos + np.array([dx, dy, dz])
            placed.add(neighbor)
            queue.append(neighbor)
    
    return coords


def generate_polymer_coords_robust(graph, target_bond_length=0.9):
    """
    Generate polymer coordinates with robust handling of various topologies.
    """
    n_atoms = graph.number_of_nodes()
    
    # Get connected components
    components = list(nx.connected_components(graph))
    
    if len(components) == 1:
        component = components[0]
        
        # Check if it's a simple cycle
        if graph.number_of_edges() == n_atoms and all(graph.degree(n) == 2 for n in graph.nodes()):
            # Perfect ring
            coords_array = place_ring(n_atoms, target_bond_length)
            coords = {i: coords_array[i] for i in range(n_atoms)}
        else:
            # Single connected component, not a simple ring
            coords = place_component_bfs(graph, component, target_bond_length)
    else:
        # Multiple components - handle each separately
        coords = {}
        offset = np.array([0.0, 0.0, 0.0])
        
        for comp_idx, component in enumerate(components):
            # Create subgraph for this component
            subgraph = graph.subgraph(component)
            
            # Check if component is a ring
            if (subgraph.number_of_edges() == len(component) and 
                all(subgraph.degree(n) == 2 for n in component)):
                # Place as ring
                comp_coords_array = place_ring(len(component), target_bond_length)
                comp_coords = {list(component)[i]: comp_coords_array[i] + offset 
                             for i in range(len(component))}
            else:
                # Place using BFS
                comp_coords = place_component_bfs(graph, component, target_bond_length)
                # Apply offset
                for node in comp_coords:
                    comp_coords[node] = comp_coords[node] + offset
            
            coords.update(comp_coords)
            
            # Update offset for next component
            if comp_coords:
                max_coord = max(np.linalg.norm(pos) for pos in comp_coords.values())
                offset[0] += 2 * max_coord + 5.0
    
    # Convert to array
    coords_array = np.zeros((n_atoms, 3))
    for i in range(n_atoms):
        if i in coords:
            coords_array[i] = coords[i]
    
    # Center the molecule
    coords_array -= np.mean(coords_array, axis=0)
    
    # Verify and adjust bond lengths
    max_attempts = 10
    for attempt in range(max_attempts):
        all_good = True
        
        for u, v in graph.edges():
            dist = np.linalg.norm(coords_array[u] - coords_array[v])
            if dist > 1.35:  # Too close to FENE limit
                # Pull atoms closer
                direction = (coords_array[v] - coords_array[u]) / dist
                target_dist = target_bond_length
                coords_array[v] = coords_array[u] + direction * target_dist
                all_good = False
            elif dist < 0.7:  # Too close
                # Push atoms apart
                direction = (coords_array[v] - coords_array[u]) / (dist + 1e-10)
                target_dist = target_bond_length
                coords_array[v] = coords_array[u] + direction * target_dist
                all_good = False
        
        if all_good:
            break
    
    return coords_array


def run_polymer_simulation(mol_id=0, lammps_exe='lmp', quick=False, data_dir=None):
    """Run simulation for a single polymer using graph structure."""
    
    # Find data directory
    if data_dir is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(script_dir, 'data')
    
    pickle_path = os.path.join(data_dir, 'rg2.pickle')
    csv_path = os.path.join(data_dir, 'rg2_SMILES.csv')
    
    # Load data
    print(f"Loading polymer data for molecule {mol_id}...")
    with open(pickle_path, 'rb') as f:
        x, rg2, topo_desc, topo_class, poly_param, graph_array = [pickle.load(f) for _ in range(6)]
    
    df = pd.read_csv(csv_path)
    
    # Get polymer info
    topology = df.iloc[mol_id]['Topology']
    n_atoms_csv = int(df.iloc[mol_id]['Number of nodes'])
    ref_rg2_ang = df.iloc[mol_id]['Rg2']
    ref_rg2_lj = rg2[mol_id, 0]
    
    # Get the graph object
    g_orig = graph_array[mol_id]
    
    # Check if graph has non-contiguous indices
    node_list = sorted(g_orig.nodes())
    if node_list != list(range(len(node_list))):
        print(f"  WARNING: Graph has non-contiguous node indices")
        print(f"  Remapping to contiguous indices 0 to {len(node_list)-1}")
        g = fix_graph_indices(g_orig)
    else:
        g = g_orig
    
    n_atoms = g.number_of_nodes()
    n_bonds = g.number_of_edges()
    
    # Check connectivity
    n_components = nx.number_connected_components(g)
    
    print(f"Molecule {mol_id}: {topology}")
    print(f"  Atoms: {n_atoms} (CSV: {n_atoms_csv})")
    print(f"  Bonds: {n_bonds}")
    print(f"  Connected components: {n_components}")
    print(f"  Reference Rg² = {ref_rg2_ang:.3f} Å² = {ref_rg2_lj:.3f} σ²")
    
    # Generate coordinates
    coords = generate_polymer_coords_robust(g, target_bond_length=0.9)
    
    # Verify bond lengths
    print("  Checking bond lengths...")
    bond_lengths = []
    for u, v in g.edges():
        dist = np.linalg.norm(coords[u] - coords[v])
        bond_lengths.append(dist)
    
    if bond_lengths:
        print(f"  Bond lengths: min={min(bond_lengths):.3f}, max={max(bond_lengths):.3f}, "
              f"mean={np.mean(bond_lengths):.3f}")
        if max(bond_lengths) > 1.4:
            print("  WARNING: Some bonds are close to FENE limit!")
    
    # Extract bonds from graph
    bonds = [(i+1, j+1) for i, j in g.edges()]  # 1-indexed for LAMMPS
    
    # Write LAMMPS data file
    data_file = f'polymer_{mol_id}.data'
    with open(data_file, 'w') as f:
        f.write(f"LAMMPS data file for polymer {mol_id}\n\n")
        f.write(f"{n_atoms} atoms\n")
        f.write(f"{n_bonds} bonds\n\n")  # Make sure this is integer
        f.write("1 atom types\n")
        f.write("1 bond types\n\n")
        
        # Box size
        max_coord = np.max(np.abs(coords))
        box_size = max(30.0, max_coord + 10.0)
        f.write(f"{-box_size:.3f} {box_size:.3f} xlo xhi\n")
        f.write(f"{-box_size:.3f} {box_size:.3f} ylo yhi\n")
        f.write(f"{-box_size:.3f} {box_size:.3f} zlo zhi\n\n")
        
        f.write("Masses\n\n")
        f.write("1 1.0\n\n")
        
        f.write("Atoms\n\n")
        for i in range(n_atoms):
            x, y, z = coords[i]
            # atom-ID molecule-ID atom-type x y z
            f.write(f"{i+1} 1 1 {x:.6f} {y:.6f} {z:.6f}\n")
        
        if n_bonds > 0:
            f.write("\nBonds\n\n")
            for i, (a1, a2) in enumerate(bonds):
                f.write(f"{i+1} 1 {a1} {a2}\n")
    
    # Simulation parameters
    if quick:
        equil_steps = 100000
        prod_steps = 100000
    else:
        equil_steps = 10000000
        prod_steps = 10000000
    
    # Write LAMMPS input script
    input_file = f'polymer_{mol_id}.in'
    with open(input_file, 'w') as f:
        f.write(f"""# LAMMPS input for Kremer-Grest polymer
# Molecule {mol_id}: {topology}, {n_atoms} atoms, {n_bonds} bonds

# Initialization
units lj
atom_style bond
boundary p p p

# Read structure
read_data {data_file}

# Force field - WCA potential
pair_style lj/cut 1.122462048
pair_coeff 1 1 1.0 1.0 1.122462048
pair_modify shift yes

# FENE bonds (only if bonds exist)
""")
        if n_bonds > 0:
            f.write("""bond_style fene
bond_coeff 1 30.0 1.5 1.0 1.0
special_bonds fene
""")
        
        f.write(f"""
# Settings
neighbor 0.3 bin
neigh_modify every 1 delay 0 check yes

# Initial minimization
thermo 100
minimize 1.0e-6 1.0e-8 100000 1000000

# Set velocities
velocity all create 1.0 {np.random.randint(10000, 99999)} dist gaussian

# Langevin dynamics
fix 1 all langevin 1.0 1.0 10.0 {np.random.randint(10000, 99999)}
fix 2 all nve

timestep 0.001

# Compute radius of gyration
compute rg all gyration

# Need to run 0 steps before evaluating variable
run 0

# Now define variable
variable rgsq equal c_rg*c_rg

# Equilibration
thermo 10000
thermo_style custom step temp pe ke etotal press c_rg
print "Starting equilibration..."
run {equil_steps}

# Production run
reset_timestep 0
fix rg_out all ave/time 2000 1 2000 v_rgsq file rg_{mol_id}.dat
print "Starting production..."
run {prod_steps}

# Final statistics
variable rg_final equal c_rg
fix rg_final_ave all ave/time 1 100000 100000 v_rgsq ave one file rg_final_{mol_id}.dat
run 100000

print "Final Rg = ${{rg_final}} sigma"
""")
    
    # Run LAMMPS
    print(f"Running LAMMPS simulation...")
    print(f"  Equilibration: {equil_steps} steps")
    print(f"  Production: {prod_steps} steps")
    start_time = time.time()
    
    log_file = f'polymer_{mol_id}.log'
    with open(log_file, 'w') as log:
        result = subprocess.run(
            [lammps_exe, '-in', input_file],
            stdout=log,
            stderr=subprocess.STDOUT
        )
    
    runtime = time.time() - start_time
    print(f"Simulation completed in {runtime:.1f} seconds")
    
    if result.returncode == 0:
        # Extract results
        try:
            # Read production data
            data = np.loadtxt(f'rg_{mol_id}.dat', skiprows=2)
            if data.size > 0:
                if data.ndim == 1:
                    rg2_sim = data[1]
                    rg2_std = 0.0
                else:
                    # Average last half of production
                    mid = len(data) // 2
                    rg2_sim = np.mean(data[mid:, 1])
                    rg2_std = np.std(data[mid:, 1])
                
                # Read final average
                try:
                    final_data = np.loadtxt(f'rg_final_{mol_id}.dat', skiprows=2)
                    if final_data.size >= 2:
                        rg2_final = final_data[1] if final_data.ndim == 1 else final_data[-1, 1]
                    else:
                        rg2_final = rg2_sim
                except:
                    rg2_final = rg2_sim
                
                # Calculate conversion factor
                sigma = np.sqrt(ref_rg2_ang / ref_rg2_lj)
                
                # Calculate errors
                error_lj = abs(rg2_sim - ref_rg2_lj) / ref_rg2_lj * 100
                
                print(f"\nResults:")
                print(f"  Simulated Rg² = {rg2_sim:.3f} ± {rg2_std:.3f} σ²")
                print(f"  Final Rg² = {rg2_final:.3f} σ²")
                print(f"  Reference Rg² = {ref_rg2_lj:.3f} σ²")
                print(f"  Error = {error_lj:.1f}%")
                print(f"  Conversion σ = {sigma:.3f} Å")
                
                return {
                    'mol_id': mol_id,
                    'topology': topology,
                    'n_atoms': n_atoms,
                    'n_bonds': n_bonds,
                    'rg2_sim': rg2_sim,
                    'rg2_std': rg2_std,
                    'rg2_final': rg2_final,
                    'rg2_ref_lj': ref_rg2_lj,
                    'rg2_ref_ang': ref_rg2_ang,
                    'error_percent': error_lj,
                    'sigma': sigma,
                    'runtime': runtime,
                    'success': True
                }
            else:
                print("ERROR: No Rg data found")
        except Exception as e:
            print(f"ERROR extracting results: {e}")
    else:
        print(f"ERROR: LAMMPS failed with code {result.returncode}")
        with open(log_file, 'r') as f:
            lines = f.readlines()
            print("\nLast 20 lines of log:")
            for line in lines[-20:]:
                print(f"  {line.strip()}")
    
    return {'mol_id': mol_id, 'success': False}


if __name__ == '__main__':
    mol_id = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    lammps_exe = 'lmp'
    quick = False
    
    for arg in sys.argv[2:]:
        if arg == '--quick':
            quick = True
        else:
            lammps_exe = arg
    
    # Create work directory
    now = datetime.now()
    timestamp = f"{now.year}_{now.month}_{now.day}_{now.hour:02d}{now.minute:02d}"
    work_dir = os.path.join(timestamp, f'polymer{mol_id}')
    os.makedirs(work_dir, exist_ok=True)
    
    original_dir = os.getcwd()
    os.chdir(work_dir)
    
    try:
        result = run_polymer_simulation(mol_id, lammps_exe, quick)
        
        if result and result.get('success'):
            print("\nSimulation successful!")
            print(f"Results saved in: {work_dir}/")
            
            # Save CSV
            csv_data = {
                'polymer_idx': [mol_id],
                'reference_rg2': [result['rg2_ref_ang']],
                'rg2': [result['rg2_sim'] * result['sigma']**2],
                'error': [result['error_percent']]
            }
            df = pd.DataFrame(csv_data)
            csv_file = f'polymer_{mol_id}_results.csv'
            df.to_csv(csv_file, index=False)
            print(f"CSV results saved to: {csv_file}")
        else:
            print("\nSimulation failed!")
    finally:
        os.chdir(original_dir)