#!/usr/bin/env python3
"""
Generate a Newick format phylogeny topology from FST distance matrix.

This script reads a pairwise FST file where the last row contains mean FST values,
constructs a distance matrix, and generates a phylogenetic tree using UPGMA clustering.
"""

import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree
import argparse
import sys


def parse_fst_file(fst_file):
    """
    Parse the FST file and extract population names and FST values.
    
    Args:
        fst_file (str): Path to the FST file
        
    Returns:
        tuple: (populations, fst_matrix) where populations is a list of population names
               and fst_matrix is a symmetric distance matrix
    """
    # Read the file
    df = pd.read_csv(fst_file, sep='\t', index_col=0)
    
    # Get the fst_mean row
    fst_row = df.loc['fst_mean']
    
    # Parse population pairs and FST values
    populations = set()
    fst_pairs = {}
    
    for pair, fst_value in fst_row.items():
        if ':' in pair:
            pop1, pop2 = pair.split(':')
            populations.add(pop1)
            populations.add(pop2)
            fst_pairs[(pop1, pop2)] = float(fst_value)
            fst_pairs[(pop2, pop1)] = float(fst_value)  # Make symmetric
    
    # Convert to sorted list for consistent indexing
    populations = sorted(list(populations))
    n_pops = len(populations)
    
    # Create distance matrix
    fst_matrix = np.zeros((n_pops, n_pops))
    
    for i, pop1 in enumerate(populations):
        for j, pop2 in enumerate(populations):
            if i == j:
                fst_matrix[i, j] = 0.0
            else:
                if (pop1, pop2) in fst_pairs:
                    fst_matrix[i, j] = fst_pairs[(pop1, pop2)]
                else:
                    # If pair not found, use a large distance or handle missing data
                    print(f"Warning: Missing FST value for {pop1}:{pop2}")
                    fst_matrix[i, j] = np.nan
    
    return populations, fst_matrix


def build_newick_tree(populations, distance_matrix, method='average'):
    """
    Build a Newick format tree from distance matrix using hierarchical clustering.
    
    Args:
        populations (list): List of population names
        distance_matrix (numpy.ndarray): Symmetric distance matrix
        method (str): Linkage method for clustering ('average' for UPGMA)
        
    Returns:
        str: Newick format tree string
    """
    # Check for NaN values and fill with maximum distance if needed
    if np.any(np.isnan(distance_matrix)):
        max_dist = np.nanmax(distance_matrix)
        distance_matrix = np.nan_to_num(distance_matrix, nan=max_dist * 1.1)
    
    # Convert to condensed distance matrix (upper triangle)
    from scipy.spatial.distance import squareform
    condensed_distances = squareform(distance_matrix)
    
    # Perform hierarchical clustering
    linkage_matrix = linkage(condensed_distances, method=method)
    
    # Convert to tree
    tree = to_tree(linkage_matrix, rd=False)
    
    def get_newick(node, parent_dist=0.0):
        """
        Recursively build Newick string from scipy tree.
        """
        if node.is_leaf():
            return f"{populations[node.id]}:{parent_dist - node.dist:.6f}"
        else:
            left_newick = get_newick(node.get_left(), node.dist)
            right_newick = get_newick(node.get_right(), node.dist)
            return f"({left_newick},{right_newick}):{parent_dist - node.dist:.6f}"
    
    # Generate Newick string
    newick = get_newick(tree, tree.dist)
    return newick + ";"


def main():
    parser = argparse.ArgumentParser(description='Generate Newick phylogeny from FST data')
    parser.add_argument('fst_file', help='Input FST file (TSV format)')
    parser.add_argument('-o', '--output', help='Output Newick file (default: stdout)')
    parser.add_argument('-m', '--method', default='average', 
                       choices=['single', 'complete', 'average', 'weighted'],
                       help='Clustering method (default: average/UPGMA)')
    
    args = parser.parse_args()
    
    try:
        # Parse FST file
        print(f"Reading FST data from {args.fst_file}...", file=sys.stderr)
        populations, fst_matrix = parse_fst_file(args.fst_file)
        
        print(f"Found {len(populations)} populations", file=sys.stderr)
        
        # Build phylogenetic tree
        print("Building phylogenetic tree...", file=sys.stderr)
        newick_tree = build_newick_tree(populations, fst_matrix, args.method)
        
        # Output result
        if args.output:
            with open(args.output, 'w') as f:
                f.write(newick_tree + '\n')
            print(f"Newick tree written to {args.output}", file=sys.stderr)
        else:
            print(newick_tree)
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()