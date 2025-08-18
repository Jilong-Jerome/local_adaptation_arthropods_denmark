#!/usr/bin/env python3
"""
Phylogeny-based analysis of selection regions
Creates phylogeny from Fst distances and classifies selection region distributions

Usage:
    python phylogeny_selection_analysis.py <fst_file> <selection_regions_file> <output_folder> [min_length]

Arguments:
    fst_file: Path to FST file with population pairwise distances
    selection_regions_file: Path to selection regions TSV file
    output_folder: Directory for output files (created if doesn't exist)
    min_length: Minimum region length in bp (default: 1000)
"""

import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Phylo
from io import StringIO
import warnings
import sys
import os
import argparse
from pathlib import Path
warnings.filterwarnings('ignore')

def parse_fst_data(fst_file):
    """Parse Fst file to extract population names and distance matrix"""
    df = pd.read_csv(fst_file, sep='\t')
    
    print(f"FST file has {len(df)} rows")
    print(f"Available rows: {df.iloc[:, 0].tolist()}")
    
    # Find the fst_mean row dynamically
    fst_mean_row = None
    for i, row in df.iterrows():
        if row.iloc[0] == 'fst_mean':
            fst_mean_row = row
            print(f"Found 'fst_mean' row at index {i}")
            break
    
    if fst_mean_row is None:
        raise ValueError("Could not find 'fst_mean' row in FST file. Available rows: " + str(df.iloc[:, 0].tolist()))
    
    # Get population comparisons from column names (skip first column 'name')
    comparisons = fst_mean_row.index[1:]
    fst_values = fst_mean_row.iloc[1:].astype(float)
    
    print(f"Found {len(fst_values)} FST comparisons")
    
    # Extract unique population names
    populations = set()
    for comp in comparisons:
        if ':' in comp:
            pop1, pop2 = comp.split(':')
            populations.add(pop1)
            populations.add(pop2)
    
    populations = sorted(list(populations))
    n_pops = len(populations)
    
    # Create symmetric distance matrix
    dist_matrix = np.zeros((n_pops, n_pops))
    pop_to_idx = {pop: i for i, pop in enumerate(populations)}
    
    for comp, fst_val in zip(comparisons, fst_values):
        if ':' in comp and not pd.isna(fst_val):
            pop1, pop2 = comp.split(':')
            if pop1 in pop_to_idx and pop2 in pop_to_idx:
                i, j = pop_to_idx[pop1], pop_to_idx[pop2]
                dist_matrix[i, j] = fst_val
                dist_matrix[j, i] = fst_val
    
    print(f"Built distance matrix for {len(populations)} populations")
    return populations, dist_matrix

def create_phylogeny(populations, dist_matrix, method='average'):
    """Create phylogeny from distance matrix using hierarchical clustering"""
    # Convert to condensed distance matrix for scipy
    condensed_dist = squareform(dist_matrix)
    
    # Perform hierarchical clustering
    linkage_matrix = linkage(condensed_dist, method=method)
    
    return linkage_matrix

def newick_from_linkage(linkage_matrix, populations):
    """Convert linkage matrix to Newick format"""
    def build_newick(node, parent_dist=0.0):
        if node < len(populations):
            # Leaf node
            return populations[node]
        else:
            # Internal node
            cluster_idx = node - len(populations)
            left_child = int(linkage_matrix[cluster_idx, 0])
            right_child = int(linkage_matrix[cluster_idx, 1])
            distance = linkage_matrix[cluster_idx, 2]
            
            left_newick = build_newick(left_child)
            right_newick = build_newick(right_child)
            
            return f"({left_newick}:{distance/2:.6f},{right_newick}:{distance/2:.6f})"
    
    root_node = len(populations) + len(linkage_matrix) - 1
    newick_string = build_newick(root_node) + ";"
    return newick_string

def parse_selection_regions(regions_file):
    """Parse selection regions file"""
    df = pd.read_csv(regions_file, sep='\t')
    
    # Calculate region length
    df['region_length'] = df['end'] - df['start'] + 1
    
    # Group by region_id to get all populations per region
    regions_by_id = df.groupby('region_id').agg({
        'species': 'first',
        'chromosome': 'first', 
        'start': 'first',
        'end': 'first',
        'region_length': 'first',
        'population_id': lambda x: list(x),
        'num_populations': 'first'
    }).reset_index()
    
    return df, regions_by_id

def map_selection_to_fst_population(selection_pop_name):
    """Map selection region population names to Fst population names"""
    # Remove EntNic_ prefix from selection region population names
    if selection_pop_name.startswith('EntNic_'):
        return selection_pop_name.replace('EntNic_', '')
    return selection_pop_name

def get_tree_tips_from_newick(newick_string):
    """Extract tip names from Newick string in tree order"""
    tree = Phylo.read(StringIO(newick_string), "newick")
    tips = [tip.name for tip in tree.get_terminals()]
    return tips

def classify_region_distribution(region_pops, tree_tips, tree, fst_populations):
    """
    Classify selection region distribution on phylogeny
    
    Returns:
    - monophyletic: all populations form a single clade
    - almost_monophyletic: populations form 2-3 clades 
    - multiphyletic: populations scattered across >3 clades
    """
    if len(region_pops) <= 1:
        return "monophyletic", 1, len(region_pops)
    
    # Map selection region population names to Fst population names
    mapped_pops = []
    unmapped_pops = []
    for pop in region_pops:
        mapped_pop = map_selection_to_fst_population(pop)
        if mapped_pop in fst_populations:
            mapped_pops.append(mapped_pop)
        else:
            unmapped_pops.append(pop)
    
    # Debug info: report unmapped populations
    if unmapped_pops:
        print(f"Warning: Could not map {len(unmapped_pops)} populations to Fst data: {unmapped_pops[:3]}{'...' if len(unmapped_pops) > 3 else ''}")
    
    if len(mapped_pops) <= 1:
        return "monophyletic", 1, len(region_pops)
    
    # Find positions of mapped populations in tree
    pop_positions = []
    for pop in mapped_pops:
        if pop in tree_tips:
            pop_positions.append(tree_tips.index(pop))
    
    if len(pop_positions) <= 1:
        return "monophyletic", 1, len(region_pops)
    
    pop_positions.sort()
    
    # Count gaps between consecutive populations
    gaps = []
    for i in range(len(pop_positions) - 1):
        gap = pop_positions[i+1] - pop_positions[i] - 1
        gaps.append(gap)
    
    # Count number of clades (groups of consecutive populations)
    clades = 1
    for gap in gaps:
        if gap > 0:  # Gap indicates separate clade
            clades += 1
    
    if clades == 1:
        return "monophyletic", clades, len(region_pops)
    elif clades <= 3:
        return "almost_monophyletic", clades, len(region_pops)
    else:
        return "multiphyletic", clades, len(region_pops)

def validate_population_mapping(regions_df, fst_populations):
    """Validate and report population mapping between files"""
    selection_pops = set(regions_df['population_id'].unique())
    
    print(f"\nPopulation mapping validation:")
    print(f"Selection regions file: {len(selection_pops)} unique populations")
    print(f"Fst file: {len(fst_populations)} populations")
    
    # Test mapping
    mapped_count = 0
    unmapped_pops = []
    
    for sel_pop in selection_pops:
        mapped_pop = map_selection_to_fst_population(sel_pop)
        if mapped_pop in fst_populations:
            mapped_count += 1
        else:
            unmapped_pops.append((sel_pop, mapped_pop))
    
    print(f"Successfully mapped: {mapped_count}/{len(selection_pops)} populations ({mapped_count/len(selection_pops)*100:.1f}%)")
    
    if unmapped_pops:
        print(f"Unmapped populations (showing first 10):")
        for sel_pop, mapped_pop in unmapped_pops[:10]:
            print(f"  {sel_pop} -> {mapped_pop} (not found in Fst)")
        
        print(f"\nFirst 10 Fst populations for reference:")
        for pop in sorted(list(fst_populations))[:10]:
            print(f"  {pop}")
    
    return mapped_count, len(selection_pops)

def analyze_selection_phylogeny(fst_file, regions_file, min_length=1000):
    """Main analysis function"""
    print("Parsing Fst data...")
    populations, dist_matrix = parse_fst_data(fst_file)
    print(f"Found {len(populations)} populations")
    
    print("Parsing selection regions...")
    regions_df, regions_by_id = parse_selection_regions(regions_file)
    
    # Validate population mapping
    validate_population_mapping(regions_df, set(populations))
    
    print("Creating phylogeny...")
    linkage_matrix = create_phylogeny(populations, dist_matrix)
    newick_string = newick_from_linkage(linkage_matrix, populations)
    
    # Filter by minimum length
    long_regions = regions_by_id[regions_by_id['region_length'] >= min_length].copy()
    print(f"Found {len(long_regions)} regions >= {min_length}bp")
    
    print("Classifying region distributions...")
    tree = Phylo.read(StringIO(newick_string), "newick")
    tree_tips = get_tree_tips_from_newick(newick_string)
    
    classifications = []
    clade_counts = []
    total_pops = []
    
    # Print mapping summary
    print(f"Fst populations: {len(populations)}")
    print(f"Tree tips: {len(tree_tips)}")
    
    for _, region in long_regions.iterrows():
        region_pops = region['population_id']
        classification, n_clades, total_region_pops = classify_region_distribution(region_pops, tree_tips, tree, set(populations))
        classifications.append(classification)
        clade_counts.append(n_clades)
        total_pops.append(total_region_pops)
    
    long_regions['phylo_distribution'] = classifications
    long_regions['n_clades'] = clade_counts
    long_regions['total_populations'] = total_pops
    
    # Summary statistics
    print("\nDistribution classification summary:")
    class_counts = long_regions['phylo_distribution'].value_counts()
    for class_type, count in class_counts.items():
        percentage = (count / len(long_regions)) * 100
        print(f"{class_type}: {count} regions ({percentage:.1f}%)")
    
    return {
        'populations': populations,
        'distance_matrix': dist_matrix,
        'linkage_matrix': linkage_matrix, 
        'newick_string': newick_string,
        'tree': tree,
        'regions_classified': long_regions,
        'classification_summary': class_counts
    }

def plot_phylogeny_with_regions(results, output_prefix):
    """Create visualizations"""
    
    # 1. Plot phylogeny
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
    
    # Dendrogram
    dendrogram(results['linkage_matrix'], 
               labels=results['populations'],
               ax=ax1, 
               orientation='top',
               leaf_rotation=90)
    ax1.set_title('Population Phylogeny (Fst distances)', fontsize=14)
    ax1.set_xlabel('Populations')
    ax1.set_ylabel('Fst Distance')
    
    # Distance matrix heatmap
    im = ax2.imshow(results['distance_matrix'], cmap='viridis', aspect='auto')
    ax2.set_xticks(range(len(results['populations'])))
    ax2.set_yticks(range(len(results['populations'])))
    ax2.set_xticklabels(results['populations'], rotation=90, fontsize=8)
    ax2.set_yticklabels(results['populations'], fontsize=8)
    ax2.set_title('Fst Distance Matrix', fontsize=14)
    
    plt.colorbar(im, ax=ax2, label='Fst')
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_phylogeny.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Plot classification summary
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    class_counts = results['classification_summary']
    colors = ['#2E86AB', '#A23B72', '#F18F01']  # Blue, purple, orange
    
    # Pie chart
    ax1.pie(class_counts.values, labels=class_counts.index, autopct='%1.1f%%',
            colors=colors[:len(class_counts)])
    ax1.set_title('Selection Region Distribution Types')
    
    # Bar chart
    bars = ax2.bar(class_counts.index, class_counts.values, color=colors[:len(class_counts)])
    ax2.set_title('Selection Region Counts by Distribution Type')
    ax2.set_ylabel('Number of Regions')
    ax2.set_xlabel('Distribution Type')
    
    # Add count labels on bars
    for bar in bars:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{int(height)}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_classification_summary.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plots saved as {output_prefix}_phylogeny.pdf and {output_prefix}_classification_summary.pdf")

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Phylogenetic analysis of selection regions based on FST distances',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python phylogeny_selection_analysis.py data/fst.tsv data/selection_regions.tsv results/ 1000
    python phylogeny_selection_analysis.py fst.tsv regions.tsv output/
        """)
    
    parser.add_argument('fst_file', 
                       help='Path to FST file with population pairwise distances')
    parser.add_argument('selection_regions_file',
                       help='Path to selection regions TSV file')
    parser.add_argument('output_folder',
                       help='Directory for output files (created if does not exist)')
    parser.add_argument('min_length', nargs='?', type=int, default=1000,
                       help='Minimum region length in bp (default: 1000)')
    
    return parser.parse_args()

def validate_inputs(args):
    """Validate input files and create output directory"""
    # Check input files exist
    if not os.path.isfile(args.fst_file):
        raise FileNotFoundError(f"FST file not found: {args.fst_file}")
    
    if not os.path.isfile(args.selection_regions_file):
        raise FileNotFoundError(f"Selection regions file not found: {args.selection_regions_file}")
    
    # Create output directory if it doesn't exist
    output_path = Path(args.output_folder)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print(f"Input validation:")
    print(f"- FST file: {args.fst_file} ✓")
    print(f"- Selection regions file: {args.selection_regions_file} ✓")
    print(f"- Output directory: {args.output_folder} ✓")
    print(f"- Minimum region length: {args.min_length} bp\n")
    
    return str(output_path.absolute())

def main():
    """Main execution"""
    # Parse command line arguments
    args = parse_arguments()
    
    # Validate inputs and create output directory
    try:
        output_dir = validate_inputs(args)
    except (FileNotFoundError, PermissionError) as e:
        print(f"Error: {e}")
        sys.exit(1)
    
    # Create output prefix
    output_prefix = os.path.join(output_dir, "phylogeny_selection")
    
    print(f"Starting phylogenetic analysis...")
    print(f"Output files will be saved with prefix: {output_prefix}\n")
    
    try:
        # Run analysis
        results = analyze_selection_phylogeny(args.fst_file, args.selection_regions_file, 
                                            min_length=args.min_length)
        
        # Create visualizations
        plot_phylogeny_with_regions(results, output_prefix)
        
        # Save results
        results['regions_classified'].to_csv(f"{output_prefix}_classified_regions.csv", index=False)
        
        # Save tree in Newick format
        with open(f"{output_prefix}_tree.newick", 'w') as f:
            f.write(results['newick_string'])
        
        # Save population mapping for reference
        fst_pops = sorted(results['populations'])
        with open(f"{output_prefix}_population_mapping.txt", 'w') as f:
            f.write("Population mapping between files:\n")
            f.write("FST file populations:\n")
            for pop in fst_pops:
                f.write(f"  {pop}\n")
            f.write(f"\nMapping rule: Remove 'EntNic_' prefix from selection region population names\n")
        
        print(f"\n🎉 Analysis completed successfully!")
        print(f"\nResults saved in: {output_dir}")
        print(f"- Classified regions: {os.path.basename(output_prefix)}_classified_regions.csv")
        print(f"- Tree: {os.path.basename(output_prefix)}_tree.newick")
        print(f"- Population mapping: {os.path.basename(output_prefix)}_population_mapping.txt")
        print(f"- Visualizations: {os.path.basename(output_prefix)}_phylogeny.pdf")
        print(f"- Classification summary: {os.path.basename(output_prefix)}_classification_summary.pdf")
        
    except Exception as e:
        print(f"\n❌ Error during analysis: {e}")
        print(f"Please check your input files and try again.")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python phylogeny_selection_analysis.py <fst_file> <selection_regions_file> <output_folder> [min_length]")
        print("\nRequired arguments:")
        print("  fst_file: Path to FST file with population pairwise distances")
        print("  selection_regions_file: Path to selection regions TSV file")
        print("  output_folder: Directory for output files (created if doesn't exist)")
        print("\nOptional arguments:")
        print("  min_length: Minimum region length in bp (default: 1000)")
        print("\nExample:")
        print("  python phylogeny_selection_analysis.py data/fst.tsv data/regions.tsv results/ 1000")
        sys.exit(1)
    
    main()
