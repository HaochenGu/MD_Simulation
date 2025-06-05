#!/usr/bin/env python3
"""
Analyze graph node indices in the polymer dataset to understand the indexing issue.
"""

import numpy as np
import pickle
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import os

def analyze_graph_indices():
    """Analyze node indices across all molecules in the dataset."""
    
    # Load data
    print("Loading data...")
    with open('data/rg2.pickle', 'rb') as f:
        x, rg2, topo_desc, topo_class, poly_param, graph_array = [pickle.load(f) for _ in range(6)]
    
    df = pd.read_csv('data/rg2_SMILES.csv')
    
    print(f"Total molecules: {len(graph_array)}")
    print(f"X array shape: {x.shape}")
    print(f"Graph array length: {len(graph_array)}")
    
    # Analyze each molecule
    results = []
    problematic_mols = []
    
    for mol_id in range(len(graph_array)):
        g = graph_array[mol_id]
        if isinstance(g, nx.Graph):
            n_nodes = g.number_of_nodes()
            n_edges = g.number_of_edges()
            
            if n_nodes > 0:
                nodes = sorted(list(g.nodes()))
                min_idx = min(nodes)
                max_idx = max(nodes)
                
                # Check if contiguous
                expected = list(range(n_nodes))
                is_contiguous = (nodes == expected)
                
                # Check if it exceeds matrix size
                exceeds_bounds = max_idx >= x.shape[1]
                
                topology = df.iloc[mol_id]['Topology']
                
                result = {
                    'mol_id': mol_id,
                    'topology': topology,
                    'n_nodes': n_nodes,
                    'n_edges': n_edges,
                    'min_idx': min_idx,
                    'max_idx': max_idx,
                    'is_contiguous': is_contiguous,
                    'exceeds_bounds': exceeds_bounds,
                    'max_expected': n_nodes - 1
                }
                results.append(result)
                
                if not is_contiguous or exceeds_bounds:
                    problematic_mols.append(mol_id)
    
    # Convert to DataFrame
    df_results = pd.DataFrame(results)
    
    # Print summary
    print("\n" + "="*60)
    print("GRAPH INDEX ANALYSIS SUMMARY")
    print("="*60)
    
    print(f"\nTotal molecules analyzed: {len(df_results)}")
    print(f"Molecules with non-contiguous indices: {(~df_results['is_contiguous']).sum()}")
    print(f"Molecules exceeding bounds: {df_results['exceeds_bounds'].sum()}")
    print(f"Total problematic molecules: {len(problematic_mols)}")
    
    # Show ranges of problematic molecules
    if problematic_mols:
        print(f"\nProblematic molecule ranges:")
        ranges = []
        start = problematic_mols[0]
        end = start
        
        for mol_id in problematic_mols[1:]:
            if mol_id == end + 1:
                end = mol_id
            else:
                ranges.append((start, end))
                start = mol_id
                end = mol_id
        ranges.append((start, end))
        
        for start, end in ranges:
            if start == end:
                print(f"  Molecule {start}")
            else:
                print(f"  Molecules {start}-{end}")
    
    # Statistics by topology
    print("\nProblematic molecules by topology:")
    for topo in sorted(df_results['topology'].unique()):
        topo_data = df_results[df_results['topology'] == topo]
        problematic = topo_data[~topo_data['is_contiguous'] | topo_data['exceeds_bounds']]
        if len(problematic) > 0:
            print(f"  {topo}: {len(problematic)}/{len(topo_data)} molecules")
    
    # Show some examples
    print("\nExample problematic molecules:")
    examples = df_results[~df_results['is_contiguous']].head(10)
    for _, row in examples.iterrows():
        print(f"  Mol {row['mol_id']}: {row['n_nodes']} nodes, "
              f"indices {row['min_idx']}-{row['max_idx']} "
              f"(expected 0-{row['max_expected']})")
    
    # Save detailed results
    output_file = 'graph_index_analysis.csv'
    df_results.to_csv(output_file, index=False)
    print(f"\nDetailed results saved to: {output_file}")
    
    # Create visualization
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Plot 1: Maximum index vs molecule ID
    ax1.scatter(df_results['mol_id'], df_results['max_idx'], 
                c=df_results['is_contiguous'], cmap='RdYlGn', alpha=0.6)
    ax1.axhline(y=99, color='red', linestyle='--', label='Matrix size limit')
    ax1.set_xlabel('Molecule ID')
    ax1.set_ylabel('Maximum Node Index')
    ax1.set_title('Maximum Node Index by Molecule ID')
    ax1.legend()
    
    # Plot 2: Index range vs number of nodes
    ax2.scatter(df_results['n_nodes'], df_results['max_idx'] - df_results['min_idx'], 
                c=df_results['is_contiguous'], cmap='RdYlGn', alpha=0.6)
    ax2.plot([0, 100], [0, 100], 'k--', label='Expected (contiguous)')
    ax2.set_xlabel('Number of Nodes')
    ax2.set_ylabel('Index Range (max - min)')
    ax2.set_title('Index Range vs Number of Nodes')
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig('graph_index_analysis.png', dpi=150)
    print(f"Visualization saved to: graph_index_analysis.png")
    
    return df_results, problematic_mols


if __name__ == '__main__':
    df_results, problematic_mols = analyze_graph_indices()
    
    # Additional analysis: check specific problematic molecule
    if problematic_mols:
        print(f"\n\nDetailed analysis of molecule {problematic_mols[0]}:")
        
        with open('data/rg2.pickle', 'rb') as f:
            x, rg2, topo_desc, topo_class, poly_param, graph_array = [pickle.load(f) for _ in range(6)]
        
        mol_id = problematic_mols[0]
        g = graph_array[mol_id]
        
        nodes = sorted(list(g.nodes()))
        print(f"All node indices: {nodes}")
        
        # Find gaps
        gaps = []
        for i in range(1, len(nodes)):
            if nodes[i] - nodes[i-1] > 1:
                gaps.append((nodes[i-1], nodes[i]))
        
        if gaps:
            print(f"\nGaps in node indices:")
            for i, (end, start) in enumerate(gaps):
                print(f"  Gap {i+1}: {end} -> {start} (missing {start-end-1} indices)")