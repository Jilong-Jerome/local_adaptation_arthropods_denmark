#!/usr/bin/env python3
"""
Explore phylogeny-selection relationships step by step
1. Build phylogeny from Fst matrix
2. Visualize Fst matrix aligned to phylogeny 
3. Create selection region presence/absence matrix aligned to phylogeny
"""

import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
import sys
import os
import argparse
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Explore phylogeny-selection relationships',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('fst_file', help='Path to FST file')
    parser.add_argument('selection_regions_file', help='Path to selection regions TSV file')
    parser.add_argument('output_folder', help='Directory for output files')
    parser.add_argument('--min_length', type=int, default=1000, help='Minimum region length (default: 1000)')
    
    return parser.parse_args()

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
        raise ValueError("Could not find 'fst_mean' row in FST file")
    
    # Get population comparisons from column names (skip first column)
    comparisons = fst_mean_row.index[1:]
    fst_values = fst_mean_row.iloc[1:].astype(float)
    
    # Extract unique population names
    populations = set()
    for comp in comparisons:
        if ':' in comp and not pd.isna(fst_values[comp]):
            pop1, pop2 = comp.split(':')
            populations.add(pop1)
            populations.add(pop2)
    
    populations = sorted(list(populations))
    n_pops = len(populations)
    
    print(f"Found {n_pops} unique populations")
    print(f"First 10 populations: {populations[:10]}")
    
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
    
    return populations, dist_matrix, comparisons, fst_values

def build_phylogeny(populations, dist_matrix, method='average'):
    """Build phylogeny using hierarchical clustering"""
    print(f"\nBuilding phylogeny using {method} linkage...")
    
    # Convert to condensed distance matrix
    condensed_dist = squareform(dist_matrix)
    
    # Perform hierarchical clustering
    linkage_matrix = linkage(condensed_dist, method=method)
    
    # Get leaf order from dendrogram
    leaf_order = leaves_list(linkage_matrix)
    ordered_populations = [populations[i] for i in leaf_order]
    
    print(f"Phylogeny built successfully")
    print(f"Population order (first 10): {ordered_populations[:10]}")
    
    return linkage_matrix, ordered_populations, leaf_order

def parse_selection_regions(regions_file, min_length=1000):
    """Parse selection regions and create summary"""
    print(f"\nParsing selection regions (min length: {min_length}bp)...")
    
    df = pd.read_csv(regions_file, sep='\t')
    print(f"Loaded {len(df)} selection region records")
    
    # Calculate region length
    df['region_length'] = df['end'] - df['start'] + 1
    
    # Filter by minimum length
    df_filtered = df[df['region_length'] >= min_length].copy()
    print(f"After filtering: {len(df_filtered)} records >= {min_length}bp")
    
    # Group by region_id to get unique regions and their populations
    regions_summary = df_filtered.groupby('region_id').agg({
        'species': 'first',
        'chromosome': 'first',
        'start': 'first', 
        'end': 'first',
        'region_length': 'first',
        'population_id': lambda x: list(x),
        'num_populations': 'first'
    }).reset_index()
    
    print(f"Found {len(regions_summary)} unique selection regions")
    
    # Get unique populations from selection data
    all_selection_pops = set()
    for pop_list in regions_summary['population_id']:
        all_selection_pops.update(pop_list)
    
    print(f"Selection regions span {len(all_selection_pops)} populations")
    
    return df_filtered, regions_summary, sorted(list(all_selection_pops))

def map_selection_to_fst_populations(selection_pops, fst_pops):
    """Map selection region population names to FST population names"""
    print(f"\nMapping population names...")
    
    def map_name(selection_pop):
        if selection_pop.startswith('EntNic_'):
            return selection_pop.replace('EntNic_', '')
        return selection_pop
    
    mapped_count = 0
    unmapped = []
    mapping = {}
    
    for sel_pop in selection_pops:
        mapped_pop = map_name(sel_pop)
        if mapped_pop in fst_pops:
            mapping[sel_pop] = mapped_pop
            mapped_count += 1
        else:
            unmapped.append(sel_pop)
    
    print(f"Successfully mapped {mapped_count}/{len(selection_pops)} populations")
    if unmapped:
        print(f"Unmapped populations (first 5): {unmapped[:5]}")
    
    return mapping

def create_selection_presence_matrix(regions_summary, population_mapping, ordered_populations):
    """Create presence/absence matrix for selection regions"""
    print(f"\nCreating selection presence/absence matrix...")
    
    # Filter regions to only include those with mappable populations
    valid_regions = []
    for _, region in regions_summary.iterrows():
        mapped_pops = []
        for pop in region['population_id']:
            if pop in population_mapping:
                mapped_pop = population_mapping[pop]
                if mapped_pop in ordered_populations:
                    mapped_pops.append(mapped_pop)
        
        if len(mapped_pops) > 0:
            region_copy = region.copy()
            region_copy['mapped_populations'] = mapped_pops
            valid_regions.append(region_copy)
    
    print(f"Found {len(valid_regions)} regions with mappable populations")
    
    # Create presence/absence matrix
    matrix = np.zeros((len(valid_regions), len(ordered_populations)), dtype=int)
    
    region_labels = []
    for i, region in enumerate(valid_regions):
        region_id = region['region_id']
        region_labels.append(region_id)
        
        for pop in region['mapped_populations']:
            if pop in ordered_populations:
                j = ordered_populations.index(pop)
                matrix[i, j] = 1
    
    print(f"Created {matrix.shape[0]} x {matrix.shape[1]} presence/absence matrix")
    
    return matrix, region_labels, valid_regions

def plot_phylogeny_fst_alignment(linkage_matrix, ordered_populations, original_populations, 
                                 dist_matrix, output_prefix):
    """Create phylogeny + FST matrix alignment plot"""
    print(f"\nCreating phylogeny-FST alignment plot...")
    
    # Use gridspec for precise control of subplot positioning
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(2, 1, height_ratios=[1, 2], hspace=0.05)
    
    ax_dendro = fig.add_subplot(gs[0])
    ax_fst = fig.add_subplot(gs[1])
    
    # 1. Dendrogram
    dendro = dendrogram(linkage_matrix, labels=ordered_populations, ax=ax_dendro,
                       orientation='top', leaf_rotation=90)
    ax_dendro.set_title('Population Phylogeny (FST distances)', fontsize=14, fontweight='bold')
    ax_dendro.set_ylabel('FST Distance')
    ax_dendro.tick_params(axis='x', labelbottom=False)  # Hide x labels on dendrogram
    
    # Get dendrogram x-coordinates for alignment
    dendro_xlim = ax_dendro.get_xlim()
    
    # 2. Reorder FST distance matrix according to phylogeny
    print(f"Original populations order: {len(original_populations)} pops")
    print(f"Phylogeny order: {len(ordered_populations)} pops")
    
    # Create mapping from original to phylogeny order
    original_to_idx = {pop: i for i, pop in enumerate(original_populations)}
    phylo_indices = [original_to_idx[pop] for pop in ordered_populations if pop in original_to_idx]
    
    print(f"Successfully mapped {len(phylo_indices)} populations for reordering")
    
    # Reorder both rows and columns of distance matrix
    reordered_dist_matrix = dist_matrix[np.ix_(phylo_indices, phylo_indices)]
    
    # Plot FST matrix with precise extent to match dendrogram
    n_pops = len(ordered_populations)
    extent = [dendro_xlim[0], dendro_xlim[1], n_pops - 0.5, -0.5]
    
    im_fst = ax_fst.imshow(reordered_dist_matrix, cmap='viridis', aspect='auto', extent=extent, origin='lower')
    
    # Set ticks to match dendrogram positions
    dendro_positions = [dendro_xlim[0] + (dendro_xlim[1] - dendro_xlim[0]) * (i + 0.5) / n_pops 
                       for i in range(n_pops)]
    
    ax_fst.set_xticks(dendro_positions)
    ax_fst.set_yticks(range(n_pops))
    ax_fst.set_xticklabels(ordered_populations, rotation=90, fontsize=10)
    ax_fst.set_yticklabels(ordered_populations, fontsize=10)
    ax_fst.set_title('FST Distance Matrix (aligned to phylogeny)', fontsize=14, fontweight='bold')
    ax_fst.set_xlabel('Populations')
    ax_fst.set_ylabel('Populations')
    
    # Set x-axis limits to match dendrogram exactly
    ax_fst.set_xlim(dendro_xlim)
    
    # Add colorbar
    plt.colorbar(im_fst, ax=ax_fst, label='FST Distance', shrink=0.8)
    
    plt.savefig(f"{output_prefix}_phylogeny_fst.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved phylogeny-FST plot: {output_prefix}_phylogeny_fst.pdf")

def plot_phylogeny_selection_alignment(linkage_matrix, ordered_populations, 
                                      selection_matrix, region_labels, output_prefix):
    """Create phylogeny + selection regions alignment plot"""
    print(f"\nCreating phylogeny-selection alignment plot...")
    
    if len(selection_matrix) == 0:
        print("No selection regions to plot")
        return
    
    # Use gridspec for precise control of subplot positioning
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(2, 1, height_ratios=[1, 3], hspace=0.05)
    
    ax_dendro = fig.add_subplot(gs[0])
    ax_sel = fig.add_subplot(gs[1])
    
    # 1. Dendrogram
    dendro = dendrogram(linkage_matrix, labels=ordered_populations, ax=ax_dendro,
                       orientation='top', leaf_rotation=90)
    ax_dendro.set_title('Population Phylogeny (FST distances)', fontsize=14, fontweight='bold')
    ax_dendro.set_ylabel('FST Distance')
    ax_dendro.tick_params(axis='x', labelbottom=False)  # Hide x labels on dendrogram
    
    # Get dendrogram x-coordinates for alignment
    dendro_xlim = ax_dendro.get_xlim()
    
    # 2. Selection regions presence/absence matrix
    print(f"Selection matrix shape: {selection_matrix.shape}")
    print(f"Populations in matrix: {len(ordered_populations)}")
    
    n_pops = len(ordered_populations)
    n_regions = len(region_labels)
    extent = [dendro_xlim[0], dendro_xlim[1], -0.5, n_regions - 0.5]
    
    im_sel = ax_sel.imshow(selection_matrix, cmap='RdBu_r', aspect='auto', 
                          vmin=0, vmax=1, extent=extent)
    
    # Set ticks to match dendrogram positions
    dendro_positions = [dendro_xlim[0] + (dendro_xlim[1] - dendro_xlim[0]) * (i + 0.5) / n_pops 
                       for i in range(n_pops)]
    
    ax_sel.set_xticks(dendro_positions)
    ax_sel.set_xticklabels(ordered_populations, rotation=90, fontsize=10)
    ax_sel.set_xlabel('Populations (phylogeny order)')
    ax_sel.set_ylabel('Selection Regions')
    ax_sel.set_title(f'Selection Regions Presence/Absence ({len(region_labels)} regions)', 
                    fontsize=14, fontweight='bold')
    
    # Set x-axis limits to match dendrogram exactly
    ax_sel.set_xlim(dendro_xlim)
    
    # Add colorbar
    cbar = plt.colorbar(im_sel, ax=ax_sel, shrink=0.8)
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(['Absent', 'Present'])
    
    # Add region labels on y-axis
    if len(region_labels) <= 100:
        ax_sel.set_yticks(range(len(region_labels)))
        ax_sel.set_yticklabels(region_labels, fontsize=6)
    else:
        # Show every nth label if too many
        step = max(1, len(region_labels) // 50)
        y_ticks = range(0, len(region_labels), step)
        ax_sel.set_yticks(y_ticks)
        ax_sel.set_yticklabels([region_labels[i] for i in y_ticks], fontsize=6)
    
    plt.savefig(f"{output_prefix}_phylogeny_selection.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved phylogeny-selection plot: {output_prefix}_phylogeny_selection.pdf")

def save_results(ordered_populations, dist_matrix, selection_matrix, region_labels, 
                valid_regions, output_prefix):
    """Save analysis results"""
    print(f"\nSaving results...")
    
    # Save population order
    pd.DataFrame({'population': ordered_populations, 
                 'phylogeny_order': range(len(ordered_populations))}).to_csv(
        f"{output_prefix}_population_order.csv", index=False)
    
    # Save selection matrix
    if len(selection_matrix) > 0:
        selection_df = pd.DataFrame(selection_matrix, 
                                   columns=ordered_populations,
                                   index=region_labels)
        selection_df.to_csv(f"{output_prefix}_selection_matrix.csv")
        
        # Save region details
        regions_details = pd.DataFrame(valid_regions)
        regions_details.to_csv(f"{output_prefix}_region_details.csv", index=False)
    
    print(f"Results saved with prefix: {output_prefix}")

def main():
    """Main execution"""
    args = parse_arguments()
    
    # Create output directory
    output_dir = Path(args.output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = output_dir / "phylogeny_exploration"
    
    try:
        # 1. Parse FST data and build phylogeny
        populations, dist_matrix, comparisons, fst_values = parse_fst_data(args.fst_file)
        linkage_matrix, ordered_populations, leaf_order = build_phylogeny(populations, dist_matrix)
        
        # 2. Parse selection regions
        regions_df, regions_summary, selection_pops = parse_selection_regions(
            args.selection_regions_file, args.min_length)
        
        # 3. Map population names
        population_mapping = map_selection_to_fst_populations(selection_pops, populations)
        
        # 4. Create selection presence/absence matrix
        selection_matrix, region_labels, valid_regions = create_selection_presence_matrix(
            regions_summary, population_mapping, ordered_populations)
        
        # 5. Create separate visualizations
        plot_phylogeny_fst_alignment(linkage_matrix, ordered_populations, populations,
                                    dist_matrix, str(output_prefix))
        
        plot_phylogeny_selection_alignment(linkage_matrix, ordered_populations,
                                          selection_matrix, region_labels, str(output_prefix))
        
        # 6. Save results
        save_results(ordered_populations, dist_matrix, selection_matrix, 
                    region_labels, valid_regions, str(output_prefix))
        
        print(f"\n🎉 Analysis completed successfully!")
        print(f"Check output files in: {output_dir}")
        print(f"\nCreated files:")
        print(f"- {output_prefix}_phylogeny_fst.pdf")
        print(f"- {output_prefix}_phylogeny_selection.pdf") 
        print(f"- {output_prefix}_population_order.csv")
        print(f"- {output_prefix}_selection_matrix.csv")
        print(f"- {output_prefix}_region_details.csv")
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()