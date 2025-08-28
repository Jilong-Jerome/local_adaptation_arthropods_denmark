#!/usr/bin/env python3
"""
Prepare aligned inputs for phylogenetic analysis.

This script takes an FST-based Newick tree and FASTA sequences, aligns their population IDs,
and outputs a topology-only Newick file paired with a FASTA file containing sequences 
with matching IDs for use with standard phylogenetic tools like IQ-TREE.
"""

import argparse
import sys
import re
from Bio import SeqIO, Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO


def load_fst_tree(newick_file):
    """
    Load FST tree from Newick file.
    
    Args:
        newick_file (str): Path to Newick file
        
    Returns:
        Bio.Phylo tree object
    """
    try:
        with open(newick_file, 'r') as f:
            newick_string = f.read().strip()
        tree = Phylo.read(StringIO(newick_string), "newick")
        return tree
    except Exception as e:
        print(f"Error loading tree: {e}", file=sys.stderr)
        sys.exit(1)


def extract_population_ids_from_fasta(fasta_file):
    """
    Extract population IDs from FASTA file headers.
    
    Args:
        fasta_file (str): Path to FASTA file
        
    Returns:
        dict: {original_sequence_id: (population_id, SeqRecord)}
    """
    sequences = {}
    population_map = {}
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Extract population ID from sequence name
        # Expected format: EntNic_POPULATION-ID_rep001
        match = re.search(r'EntNic_([^_]+)_rep\d+', record.id)
        if match:
            pop_id = match.group(1)
            sequences[record.id] = (pop_id, record)
            population_map[pop_id] = record.id
        else:
            print(f"Warning: Could not extract population ID from {record.id}", file=sys.stderr)
    
    return sequences, population_map


def get_tree_terminal_names(tree):
    """
    Get terminal (leaf) names from tree.
    
    Args:
        tree: Bio.Phylo tree object
        
    Returns:
        list: Terminal node names
    """
    terminals = tree.get_terminals()
    return [term.name for term in terminals if term.name]


def create_name_mapping(tree_names, fasta_pop_ids, prefix="T"):
    """
    Create mapping between tree names, FASTA population IDs, and standardized names.
    
    Args:
        tree_names (list): Terminal names from tree
        fasta_pop_ids (list): Population IDs from FASTA
        prefix (str): Prefix for standardized names
        
    Returns:
        tuple: (name_mapping_dict, matched_populations)
    """
    # Find intersection of tree names and FASTA population IDs
    tree_set = set(tree_names)
    fasta_set = set(fasta_pop_ids)
    matched_pops = tree_set & fasta_set
    
    print(f"Tree terminals: {len(tree_names)}", file=sys.stderr)
    print(f"FASTA populations: {len(fasta_pop_ids)}", file=sys.stderr)
    print(f"Matched populations: {len(matched_pops)}", file=sys.stderr)
    
    if len(matched_pops) == 0:
        print("Error: No matching population IDs between tree and FASTA!", file=sys.stderr)
        print(f"Tree names sample: {tree_names[:5]}", file=sys.stderr)
        print(f"FASTA IDs sample: {list(fasta_pop_ids)[:5]}", file=sys.stderr)
        sys.exit(1)
    
    # Create standardized names
    name_mapping = {}
    sorted_matched = sorted(list(matched_pops))
    
    for i, pop_id in enumerate(sorted_matched):
        # Create standardized name: T001, T002, etc.
        std_name = f"{prefix}{i+1:03d}"
        name_mapping[pop_id] = {
            'standard_name': std_name,
            'original_pop_id': pop_id,
            'in_tree': pop_id in tree_set,
            'in_fasta': pop_id in fasta_set
        }
    
    return name_mapping, matched_pops


def create_topology_tree(tree, name_mapping, matched_pops):
    """
    Create topology-only tree with standardized names.
    
    Args:
        tree: Original Bio.Phylo tree
        name_mapping (dict): Name mapping dictionary
        matched_pops (set): Set of matched population IDs
        
    Returns:
        Bio.Phylo tree object with topology only
    """
    # Create a copy of the tree using StringIO
    from io import StringIO
    tree_string = StringIO()
    Phylo.write(tree, tree_string, "newick")
    tree_string.seek(0)
    topo_tree = Phylo.read(tree_string, "newick")
    
    # Remove branch lengths and rename terminals
    terminals_to_remove = []
    
    for terminal in topo_tree.get_terminals():
        if terminal.name in matched_pops:
            # Rename to standardized name
            terminal.name = name_mapping[terminal.name]['standard_name']
            # Remove branch length
            terminal.branch_length = None
        else:
            # Mark for removal
            terminals_to_remove.append(terminal.name)
    
    # Remove unmatched terminals
    if terminals_to_remove:
        print(f"Removing {len(terminals_to_remove)} unmatched terminals from tree", file=sys.stderr)
        for term_name in terminals_to_remove:
            try:
                # Find and remove the terminal
                for terminal in topo_tree.get_terminals():
                    if terminal.name == term_name:
                        topo_tree.prune(terminal)
                        break
            except Exception as e:
                print(f"Warning: Could not remove terminal {term_name}: {e}", file=sys.stderr)
    
    # Remove all branch lengths from internal nodes
    for clade in topo_tree.find_clades():
        clade.branch_length = None
        clade.confidence = None  # Remove support values if present
    
    return topo_tree


def create_aligned_fasta(sequences, population_map, name_mapping, matched_pops):
    """
    Create FASTA file with aligned sequence IDs.
    
    Args:
        sequences (dict): Original sequences {seq_id: (pop_id, SeqRecord)}
        population_map (dict): {pop_id: original_seq_id}
        name_mapping (dict): Name mapping dictionary
        matched_pops (set): Set of matched population IDs
        
    Returns:
        list: List of SeqRecord objects with standardized names
    """
    aligned_sequences = []
    
    for pop_id in sorted(list(matched_pops)):
        if pop_id in population_map and pop_id in name_mapping:
            original_seq_id = population_map[pop_id]
            if original_seq_id in sequences:
                _, original_record = sequences[original_seq_id]
                std_name = name_mapping[pop_id]['standard_name']
                
                # Create new SeqRecord with standardized ID
                new_record = SeqRecord(
                    seq=original_record.seq,
                    id=std_name,
                    description=f"Population: {pop_id} | Original: {original_record.id}"
                )
                aligned_sequences.append(new_record)
    
    return aligned_sequences


def write_name_mapping_file(name_mapping, output_file):
    """
    Write name mapping to file for reference.
    
    Args:
        name_mapping (dict): Name mapping dictionary
        output_file (str): Output file path
    """
    with open(output_file, 'w') as f:
        f.write("Standard_Name\tOriginal_Population_ID\tIn_Tree\tIn_FASTA\n")
        for pop_id in sorted(name_mapping.keys()):
            info = name_mapping[pop_id]
            f.write(f"{info['standard_name']}\t{pop_id}\t{info['in_tree']}\t{info['in_fasta']}\n")


def main():
    parser = argparse.ArgumentParser(description='Prepare aligned phylogenetic inputs from FST tree and FASTA sequences')
    parser.add_argument('fst_newick', help='Input FST Newick tree file')
    parser.add_argument('fasta_sequences', help='Input FASTA sequence file')
    parser.add_argument('-p', '--prefix', default='phylo_aligned', 
                       help='Output file prefix (default: phylo_aligned)')
    parser.add_argument('--name-prefix', default='T', 
                       help='Prefix for standardized taxon names (default: T)')
    parser.add_argument('--keep-original-names', action='store_true',
                       help='Keep original population names instead of standardized names')
    
    args = parser.parse_args()
    
    try:
        # Load FST tree
        print("Loading FST tree...", file=sys.stderr)
        tree = load_fst_tree(args.fst_newick)
        
        # Extract population IDs from FASTA
        print("Processing FASTA sequences...", file=sys.stderr)
        sequences, population_map = extract_population_ids_from_fasta(args.fasta_sequences)
        
        # Get tree terminal names
        tree_names = get_tree_terminal_names(tree)
        fasta_pop_ids = list(population_map.keys())
        
        # Create name mapping
        print("Creating name mapping...", file=sys.stderr)
        if args.keep_original_names:
            # Use original names, just filter to matched populations
            tree_set = set(tree_names)
            fasta_set = set(fasta_pop_ids)
            matched_pops = tree_set & fasta_set
            
            name_mapping = {}
            for pop_id in matched_pops:
                name_mapping[pop_id] = {
                    'standard_name': pop_id,
                    'original_pop_id': pop_id,
                    'in_tree': True,
                    'in_fasta': True
                }
        else:
            name_mapping, matched_pops = create_name_mapping(tree_names, fasta_pop_ids, args.name_prefix)
        
        # Create topology-only tree
        print("Creating topology tree...", file=sys.stderr)
        topo_tree = create_topology_tree(tree, name_mapping, matched_pops)
        
        # Create aligned FASTA
        print("Creating aligned FASTA...", file=sys.stderr)
        aligned_sequences = create_aligned_fasta(sequences, population_map, name_mapping, matched_pops)
        
        # Write outputs
        tree_output = f"{args.prefix}_topology.newick"
        fasta_output = f"{args.prefix}_sequences.fasta"
        mapping_output = f"{args.prefix}_name_mapping.txt"
        
        # Write topology tree
        Phylo.write(topo_tree, tree_output, "newick")
        print(f"Topology tree written to: {tree_output}")
        
        # Write aligned FASTA
        SeqIO.write(aligned_sequences, fasta_output, "fasta")
        print(f"Aligned sequences written to: {fasta_output}")
        
        # Write name mapping
        write_name_mapping_file(name_mapping, mapping_output)
        print(f"Name mapping written to: {mapping_output}")
        
        print(f"\nSummary:", file=sys.stderr)
        print(f"  Matched populations: {len(matched_pops)}", file=sys.stderr)
        print(f"  Output files: {tree_output}, {fasta_output}, {mapping_output}", file=sys.stderr)
        print(f"\nReady for use with IQ-TREE:", file=sys.stderr)
        print(f"  iqtree -s {fasta_output} -g {tree_output} [other options]", file=sys.stderr)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()