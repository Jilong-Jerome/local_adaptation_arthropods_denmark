#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
import argparse
from collections import defaultdict
from scipy import stats

def parse_population_labels(label_file):
    """Parse population_label.txt to extract population groupings"""
    population_groups = {}
    
    with open(label_file, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].rstrip('\n')
        
        # Skip empty lines
        if not line.strip():
            i += 1
            continue
            
        # Handle Outgroup_1(South-East)
        if "Outgroup_1(South-East)" in line:
            i += 1
            if i < len(lines) and "BRS,VAJ,JEJ" in lines[i]:
                for pop in ["BRS", "VAJ", "JEJ"]:
                    population_groups[pop] = "Outgroup_1"
            i += 1
            continue
            
        # Main regions
        if line.strip() in ["North-West", "South-East"]:
            main_region = line.strip()
            i += 1
            
            # Process subregions within this main region
            while i < len(lines):
                subline = lines[i].rstrip('\n')
                
                # Break if we hit another main region or end
                if subline.strip() in ["North-West", "South-East"] or not subline.strip():
                    break
                    
                # Subregion name (single tab)
                if subline.startswith('\t') and not subline.startswith('\t\t'):
                    subregion = subline.strip()
                    i += 1
                    
                    # Population codes (double tab)
                    if i < len(lines) and lines[i].startswith('\t\t'):
                        pop_line = lines[i].strip()
                        populations = [p.strip() for p in pop_line.replace(',', ' ').split() if p.strip()]
                        
                        for pop in populations:
                            population_groups[pop] = f"{main_region}_{subregion}"
                        i += 1
                    continue
                    
                # Special cases within regions
                elif "Outgroup_2(East Jylland)" in subline:
                    i += 1
                    if i < len(lines) and "JHJ" in lines[i]:
                        population_groups["JHJ"] = "Outgroup_2"
                    i += 1
                    continue
                    
                i += 1
        else:
            i += 1
    
    return population_groups

def analyze_selection_regions(selection_file, population_labels, min_length=5000):
    """Analyze selection regions and create presence/absence matrix"""
    
    # Read selection regions data
    df = pd.read_csv(selection_file, sep='\t')
    
    # Filter by minimum length
    df['length'] = df['end'] - df['start'] + 1
    df_filtered = df[df['length'] >= min_length].copy()
    
    print(f"Total regions: {len(df['region_id'].unique())}")
    print(f"Regions >= {min_length}bp: {len(df_filtered['region_id'].unique())}")
    
    # Extract population from population_id (remove EntNic_ prefix and -C suffix)
    df_filtered['population'] = df_filtered['population_id'].str.replace(r'EntNic_', '', regex=True).str.replace(r'-C\d+', '', regex=True)
    
    # Map populations to labels
    df_filtered['population_label'] = df_filtered['population'].map(population_labels)
    
    # Create maximum overlap fraction matrix
    unique_regions = df_filtered['region_id'].unique()
    unique_labels = sorted(list(set(population_labels.values())))
    
    results = []
    
    for region in unique_regions:
        region_data = df_filtered[df_filtered['region_id'] == region]
        region_info = region_data.iloc[0]
        
        # Create result row with basic info
        result_row = {
            'region_id': region,
            'chromosome': region_info['chromosome'],
            'start': region_info['start'],
            'end': region_info['end'],
            'length': region_info['length'],
            'num_populations_total': len(region_data['population'].unique())
        }
        
        # Calculate maximum overlap fraction for each population label
        # Handle multiple segments per population by summing their lengths
        for label in unique_labels:
            # Get all populations in this region that belong to this label
            label_data = region_data[region_data['population_label'] == label]
            
            if len(label_data) > 0:
                # Group by population and sum segment lengths for each population
                population_totals = {}
                region_length = region_info['end'] - region_info['start'] + 1
                
                for _, row in label_data.iterrows():
                    pop = row['population']
                    segment_length = row['end'] - row['start'] + 1
                    
                    if pop not in population_totals:
                        population_totals[pop] = 0
                    population_totals[pop] += segment_length
                
                # Calculate overlap fraction for each population and take maximum
                max_fraction = 0
                for pop, total_length in population_totals.items():
                    overlap_fraction = total_length / region_length
                    max_fraction = max(max_fraction, overlap_fraction)
                
                result_row[label] = round(max_fraction, 3)
            else:
                # No populations in this label have this selection region
                result_row[label] = 0.0
            
        results.append(result_row)
    
    return pd.DataFrame(results), df_filtered

def create_long_format(results_df, population_labels):
    """Convert wide format results to long format"""
    
    # Get the population label columns (exclude the basic info columns)
    basic_columns = ['region_id', 'chromosome', 'start', 'end', 'length', 'num_populations_total']
    label_columns = [col for col in results_df.columns if col not in basic_columns]
    
    long_format_data = []
    
    for _, row in results_df.iterrows():
        for label in label_columns:
            long_row = {
                'region_id': row['region_id'],
                'chromosome': row['chromosome'],
                'start': row['start'],
                'end': row['end'],
                'length': row['length'],
                'num_populations_total': row['num_populations_total'],
                'population_label': label,
                'max_overlap_fraction': row[label]
            }
            long_format_data.append(long_row)
    
    return pd.DataFrame(long_format_data)

def categorize_presence_absence(results_df, threshold_present=0.5, threshold_absent=0.1):
    """Categorize fine-grained fractions into presence/absence/weak presence categories"""
    
    # Get the population label columns (exclude the basic info columns)
    basic_columns = ['region_id', 'chromosome', 'start', 'end', 'length', 'num_populations_total']
    label_columns = [col for col in results_df.columns if col not in basic_columns]
    
    # Create a copy for categorized results
    categorized_df = results_df.copy()
    
    # Apply categorization to each population label column
    for label in label_columns:
        categorized_df[f'{label}_category'] = categorized_df[label].apply(
            lambda x: 'present' if x >= threshold_present 
                     else 'absent' if x < threshold_absent 
                     else 'weak_present'
        )
    
    return categorized_df

def summarize_selection_distribution(categorized_df, population_labels, original_df, overlap_threshold=0.5, include_weak=True):
    """Summarize how selection regions are distributed across population labels using categories"""
    
    # Get the population label columns (exclude the basic info columns)
    basic_columns = ['region_id', 'chromosome', 'start', 'end', 'length', 'num_populations_total']
    label_columns = [col for col in categorized_df.columns if col not in basic_columns and not col.endswith('_category')]
    category_columns = [f'{label}_category' for label in label_columns]
    
    # Initialize summary data
    summary = {
        'total_regions': len(categorized_df),
        'include_weak_present': include_weak,
        'present_per_label': {},
        'weak_present_per_label': {},
        'north_west_only': {'count': 0, 'regions': []},
        'south_east_only': {'count': 0, 'regions': []},
        'outgroup_only': {'count': 0, 'regions': []},
        'shared_north_south': {'count': 0, 'regions': []},
        'shared_with_outgroups': {'count': 0, 'regions': []},
        'universal': {'count': 0, 'regions': []},
        'detailed_patterns': {},
        'length_stats': {},
        'peak_length_stats': {},
        'regions_per_population_stats': {}
    }
    
    # Categorize labels by main region
    north_west_labels = [label for label in label_columns if label.startswith('North-West')]
    south_east_labels = [label for label in label_columns if label.startswith('South-East')]
    outgroup_labels = [label for label in label_columns if label.startswith('Outgroup')]
    
    print(f"North-West labels: {north_west_labels}")
    print(f"South-East labels: {south_east_labels}")
    print(f"Outgroup labels: {outgroup_labels}")
    
    for _, row in categorized_df.iterrows():
        # Find which labels have this region as 'present' or 'weak_present'
        present_labels = [label for label in label_columns if row[f'{label}_category'] == 'present']
        weak_present_labels = [label for label in label_columns if row[f'{label}_category'] == 'weak_present']
        
        # Determine which labels to use for analysis
        if include_weak:
            analysis_labels = present_labels + weak_present_labels
        else:
            analysis_labels = present_labels
        
        # Count present and weak present regions per label
        for label in present_labels:
            if label not in summary['present_per_label']:
                summary['present_per_label'][label] = 0
            summary['present_per_label'][label] += 1
            
        for label in weak_present_labels:
            if label not in summary['weak_present_per_label']:
                summary['weak_present_per_label'][label] = 0
            summary['weak_present_per_label'][label] += 1
        
        # Categorize by distribution pattern
        present_north_west = [label for label in analysis_labels if label in north_west_labels]
        present_south_east = [label for label in analysis_labels if label in south_east_labels]
        present_outgroups = [label for label in analysis_labels if label in outgroup_labels]
        
        # Create pattern key for detailed tracking
        pattern_key = f"NW:{len(present_north_west)}_SE:{len(present_south_east)}_OUT:{len(present_outgroups)}"
        if pattern_key not in summary['detailed_patterns']:
            summary['detailed_patterns'][pattern_key] = 0
        summary['detailed_patterns'][pattern_key] += 1
        
        # Get peak length information from original data for this region
        region_original_data = original_df[original_df['region_id'] == row['region_id']]
        if include_weak:
            # Include populations that have present or weak_present for this region
            relevant_pops = [label for label in analysis_labels]
        else:
            # Only include populations that have present for this region  
            relevant_pops = present_labels
            
        # Get peak lengths for populations in this region that match our analysis criteria
        # Handle multiple segments per population by summing their lengths
        population_peak_lengths = {}  # pop -> total peak length
        
        for _, orig_row in region_original_data.iterrows():
            pop = orig_row['population_id'].replace('EntNic_', '').split('-C')[0]
            pop_label = population_labels.get(pop, 'Unknown')
            if pop_label in relevant_pops:
                segment_length = orig_row['peak_end'] - orig_row['peak_start'] + 1
                if pop not in population_peak_lengths:
                    population_peak_lengths[pop] = 0
                population_peak_lengths[pop] += segment_length
        
        # Convert to list of total peak lengths per population
        peak_lengths = list(population_peak_lengths.values())
        
        # Store region info for length analysis
        region_info = {
            'region_id': row['region_id'],
            'length': row['length'],
            'chromosome': row['chromosome'],
            'start': row['start'],
            'end': row['end'],
            'peak_lengths': peak_lengths
        }
        
        # Categorize into main distribution types and store region info
        if len(analysis_labels) == len(label_columns):
            summary['universal']['count'] += 1
            summary['universal']['regions'].append(region_info)
        elif present_outgroups and (present_north_west or present_south_east):
            summary['shared_with_outgroups']['count'] += 1
            summary['shared_with_outgroups']['regions'].append(region_info)
        elif present_north_west and present_south_east:
            summary['shared_north_south']['count'] += 1
            summary['shared_north_south']['regions'].append(region_info)
        elif present_north_west and not present_south_east and not present_outgroups:
            summary['north_west_only']['count'] += 1
            summary['north_west_only']['regions'].append(region_info)
        elif present_south_east and not present_north_west and not present_outgroups:
            summary['south_east_only']['count'] += 1
            summary['south_east_only']['regions'].append(region_info)
        elif present_outgroups and not present_north_west and not present_south_east:
            summary['outgroup_only']['count'] += 1
            summary['outgroup_only']['regions'].append(region_info)
    
    # Calculate length statistics for each category
    
    categories = ['north_west_only', 'south_east_only', 'outgroup_only', 
                 'shared_north_south', 'shared_with_outgroups', 'universal']
    
    for category in categories:
        if summary[category]['count'] > 0:
            # Region length statistics
            lengths = [r['length'] for r in summary[category]['regions']]
            summary['length_stats'][category] = {
                'count': len(lengths),
                'mean': np.mean(lengths),
                'median': np.median(lengths),
                'min': np.min(lengths),
                'max': np.max(lengths),
                'std': np.std(lengths)
            }
            
            # Peak length statistics
            all_peak_lengths = []
            for r in summary[category]['regions']:
                all_peak_lengths.extend(r['peak_lengths'])
            
            if all_peak_lengths:
                summary['peak_length_stats'][category] = {
                    'count': len(all_peak_lengths),
                    'mean': np.mean(all_peak_lengths),
                    'median': np.median(all_peak_lengths),
                    'min': np.min(all_peak_lengths),
                    'max': np.max(all_peak_lengths),
                    'std': np.std(all_peak_lengths)
                }
            else:
                summary['peak_length_stats'][category] = {
                    'count': 0, 'mean': 0, 'median': 0, 'min': 0, 'max': 0, 'std': 0
                }
        else:
            summary['length_stats'][category] = {
                'count': 0, 'mean': 0, 'median': 0, 'min': 0, 'max': 0, 'std': 0
            }
            summary['peak_length_stats'][category] = {
                'count': 0, 'mean': 0, 'median': 0, 'min': 0, 'max': 0, 'std': 0
            }
    
    # Calculate regions per population statistics for each label
    # Count how many regions each individual population has (overlap_fraction >= threshold)
    # Handle multiple segments per population by summing their lengths first
    population_region_counts = {}  # population -> count of regions
    pop_region_segments = {}  # (pop, region_id) -> total_segment_length
    pop_region_lengths = {}  # (pop, region_id) -> region_length
    
    # First pass: sum segments for each population-region combination
    for _, row in original_df.iterrows():
        pop = row['population_id'].replace('EntNic_', '').split('-C')[0]
        region_id = row['region_id']
        key = (pop, region_id)
        
        segment_length = row['end'] - row['start'] + 1
        # Get the actual region length from the full region boundaries
        # Find the region info from the categorized_df
        region_info = categorized_df[categorized_df['region_id'] == region_id].iloc[0]
        region_length = region_info['length']
        
        if key not in pop_region_segments:
            pop_region_segments[key] = 0
            pop_region_lengths[key] = region_length
        pop_region_segments[key] += segment_length
    
    # Second pass: count regions where total overlap fraction >= threshold
    for (pop, region_id), total_segment_length in pop_region_segments.items():
        region_length = pop_region_lengths[(pop, region_id)]
        overlap_fraction = total_segment_length / region_length
        
        if overlap_fraction >= overlap_threshold:
            if pop not in population_region_counts:
                population_region_counts[pop] = set()
            population_region_counts[pop].add(region_id)
    
    # Convert sets to counts
    for pop in population_region_counts:
        population_region_counts[pop] = len(population_region_counts[pop])
    
    # Calculate statistics per population label
    for label in label_columns:
        # Get all populations in this label
        populations_in_label = [pop for pop, pop_label in population_labels.items() if pop_label == label]
        
        # Get region counts for populations in this label
        region_counts = []
        for pop in populations_in_label:
            count = population_region_counts.get(pop, 0)
            region_counts.append(count)
        
        if region_counts:
            summary['regions_per_population_stats'][label] = {
                'num_populations': len(region_counts),
                'total_regions': sum(region_counts),
                'mean_regions_per_pop': np.mean(region_counts),
                'median_regions_per_pop': np.median(region_counts),
                'min_regions_per_pop': np.min(region_counts),
                'max_regions_per_pop': np.max(region_counts),
                'std_regions_per_pop': np.std(region_counts)
            }
        else:
            summary['regions_per_population_stats'][label] = {
                'num_populations': 0,
                'total_regions': 0,
                'mean_regions_per_pop': 0,
                'median_regions_per_pop': 0,
                'min_regions_per_pop': 0,
                'max_regions_per_pop': 0,
                'std_regions_per_pop': 0
            }
    
    # Perform Mann-Whitney U test between North-West and South-East populations (all regions)
    north_west_counts = []
    south_east_counts = []
    
    for label in label_columns:
        populations_in_label = [pop for pop, pop_label in population_labels.items() if pop_label == label]
        region_counts = [population_region_counts.get(pop, 0) for pop in populations_in_label]
        
        if label.startswith('North-West'):
            north_west_counts.extend(region_counts)
        elif label.startswith('South-East'):
            south_east_counts.extend(region_counts)
    
    # Perform the test if both groups have data
    if len(north_west_counts) > 0 and len(south_east_counts) > 0:
        statistic, p_value = stats.mannwhitneyu(north_west_counts, south_east_counts, alternative='two-sided')
        summary['mann_whitney_test'] = {
            'north_west_n': len(north_west_counts),
            'south_east_n': len(south_east_counts),
            'north_west_mean': np.mean(north_west_counts),
            'south_east_mean': np.mean(south_east_counts),
            'north_west_median': np.median(north_west_counts),
            'south_east_median': np.median(south_east_counts),
            'u_statistic': statistic,
            'p_value': p_value
        }
    else:
        summary['mann_whitney_test'] = None
    
    # Second Mann-Whitney U test for region-specific selections only
    # Identify North-West only and South-East only regions
    north_west_only_regions = set([r['region_id'] for r in summary['north_west_only']['regions']])
    south_east_only_regions = set([r['region_id'] for r in summary['south_east_only']['regions']])
    
    # Count region-specific selections per population
    population_specific_region_counts = {}  # population -> count of region-specific regions
    
    for _, row in original_df.iterrows():
        pop = row['population_id'].replace('EntNic_', '').split('-C')[0]
        pop_label = population_labels.get(pop, 'Unknown')
        region_id = row['region_id']
        
        # Only count if overlap_fraction >= threshold and region is region-specific
        if row['overlap_fraction'] >= overlap_threshold:
            if pop not in population_specific_region_counts:
                population_specific_region_counts[pop] = set()
            
            # Add region if it's specific to the population's geographic group
            if pop_label.startswith('North-West') and region_id in north_west_only_regions:
                population_specific_region_counts[pop].add(region_id)
            elif pop_label.startswith('South-East') and region_id in south_east_only_regions:
                population_specific_region_counts[pop].add(region_id)
    
    # Convert sets to counts
    for pop in population_specific_region_counts:
        population_specific_region_counts[pop] = len(population_specific_region_counts[pop])
    
    # Collect counts by geographic region for region-specific test
    north_west_specific_counts = []
    south_east_specific_counts = []
    
    for label in label_columns:
        populations_in_label = [pop for pop, pop_label in population_labels.items() if pop_label == label]
        
        if label.startswith('North-West'):
            for pop in populations_in_label:
                count = population_specific_region_counts.get(pop, 0)
                north_west_specific_counts.append(count)
        elif label.startswith('South-East'):
            for pop in populations_in_label:
                count = population_specific_region_counts.get(pop, 0)
                south_east_specific_counts.append(count)
    
    # Perform the region-specific test if both groups have data
    if len(north_west_specific_counts) > 0 and len(south_east_specific_counts) > 0:
        statistic_specific, p_value_specific = stats.mannwhitneyu(
            north_west_specific_counts, south_east_specific_counts, alternative='two-sided'
        )
        summary['mann_whitney_test_specific'] = {
            'north_west_n': len(north_west_specific_counts),
            'south_east_n': len(south_east_specific_counts),
            'north_west_mean': np.mean(north_west_specific_counts),
            'south_east_mean': np.mean(south_east_specific_counts),
            'north_west_median': np.median(north_west_specific_counts),
            'south_east_median': np.median(south_east_specific_counts),
            'u_statistic': statistic_specific,
            'p_value': p_value_specific,
            'north_west_only_regions': len(north_west_only_regions),
            'south_east_only_regions': len(south_east_only_regions)
        }
    else:
        summary['mann_whitney_test_specific'] = None
    
    # Third Mann-Whitney U test for peak lengths: shared vs unique regions
    # Collect peak lengths for different region types
    shared_nw_se_regions = set([r['region_id'] for r in summary['shared_north_south']['regions']])
    unique_regions = north_west_only_regions.union(south_east_only_regions)
    
    shared_peak_lengths = []
    unique_peak_lengths = []
    
    # Group data by population and region to sum segments
    population_region_peaks = {}  # (pop, region_id) -> total_peak_length
    population_region_segments = {}  # (pop, region_id) -> total_segment_length
    region_lengths = {}  # region_id -> region_length
    
    for _, row in original_df.iterrows():
        pop = row['population_id'].replace('EntNic_', '').split('-C')[0]
        pop_label = population_labels.get(pop, 'Unknown')
        region_id = row['region_id']
        key = (pop, region_id)
        
        # Store region length for overlap calculation
        if region_id not in region_lengths:
            region_lengths[region_id] = row['end'] - row['start'] + 1
        
        # Sum peak and segment lengths for each population-region combination
        peak_length = row['peak_end'] - row['peak_start'] + 1
        segment_length = row['end'] - row['start'] + 1
        
        if key not in population_region_peaks:
            population_region_peaks[key] = 0
            population_region_segments[key] = 0
        population_region_peaks[key] += peak_length
        population_region_segments[key] += segment_length
    
    # Now collect peak lengths based on summed segments
    for (pop, region_id), total_peak_length in population_region_peaks.items():
        pop_label = population_labels.get(pop, 'Unknown')
        total_segment_length = population_region_segments[(pop, region_id)]
        region_length = region_lengths[region_id]
        overlap_fraction = total_segment_length / region_length
        
        # Only include peaks from populations that meet our threshold criteria
        if overlap_fraction >= overlap_threshold:
            # Collect peak lengths based on region type
            if region_id in shared_nw_se_regions:
                # Only include peaks from NW or SE populations (exclude outgroups)
                if pop_label.startswith('North-West') or pop_label.startswith('South-East'):
                    shared_peak_lengths.append(total_peak_length)
            elif region_id in unique_regions:
                # Include peaks from any population in unique regions
                if pop_label.startswith('North-West') or pop_label.startswith('South-East'):
                    unique_peak_lengths.append(total_peak_length)
    
    # Perform the peak length test if both groups have data
    if len(shared_peak_lengths) > 0 and len(unique_peak_lengths) > 0:
        statistic_peaks, p_value_peaks = stats.mannwhitneyu(
            shared_peak_lengths, unique_peak_lengths, alternative='two-sided'
        )
        summary['mann_whitney_test_peaks'] = {
            'shared_n': len(shared_peak_lengths),
            'unique_n': len(unique_peak_lengths),
            'shared_mean': np.mean(shared_peak_lengths),
            'unique_mean': np.mean(unique_peak_lengths),
            'shared_median': np.median(shared_peak_lengths),
            'unique_median': np.median(unique_peak_lengths),
            'u_statistic': statistic_peaks,
            'p_value': p_value_peaks,
            'shared_regions_count': len(shared_nw_se_regions),
            'unique_regions_count': len(unique_regions)
        }
    else:
        summary['mann_whitney_test_peaks'] = None
    
    return summary

def write_distribution_summary(summary, output_file, overlap_threshold=0.5):
    """Write a formatted summary of selection region distribution to file"""
    
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("SELECTION REGION DISTRIBUTION SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        
        analysis_type = "present + weakly present" if summary['include_weak_present'] else "present only"
        f.write(f"Analysis based on: {analysis_type} regions\n")
        f.write(f"Overlap fraction threshold for presence: {overlap_threshold}\n")
        f.write(f"Total selection regions analyzed: {summary['total_regions']}\n\n")
        
        f.write("--- Distribution by main geographic regions ---\n")
        f.write(f"North-West only:           {summary['north_west_only']['count']:4d} ({summary['north_west_only']['count']/summary['total_regions']*100:.1f}%)\n")
        f.write(f"South-East only:           {summary['south_east_only']['count']:4d} ({summary['south_east_only']['count']/summary['total_regions']*100:.1f}%)\n")
        f.write(f"Outgroup only:             {summary['outgroup_only']['count']:4d} ({summary['outgroup_only']['count']/summary['total_regions']*100:.1f}%)\n")
        f.write(f"Shared North-West & South-East: {summary['shared_north_south']['count']:4d} ({summary['shared_north_south']['count']/summary['total_regions']*100:.1f}%)\n")
        f.write(f"Shared with outgroups:     {summary['shared_with_outgroups']['count']:4d} ({summary['shared_with_outgroups']['count']/summary['total_regions']*100:.1f}%)\n")
        f.write(f"Universal (all labels):    {summary['universal']['count']:4d} ({summary['universal']['count']/summary['total_regions']*100:.1f}%)\n\n")
        
        f.write("--- Selection region length distribution by geographic pattern ---\n")
        f.write("Category                   Count    Mean(bp)   Median(bp)   Min(bp)    Max(bp)    Std(bp)\n")
        f.write("-" * 85 + "\n")
        
        categories = [
            ('North-West only', 'north_west_only'),
            ('South-East only', 'south_east_only'), 
            ('Outgroup only', 'outgroup_only'),
            ('Shared NW & SE', 'shared_north_south'),
            ('Shared w/ outgroups', 'shared_with_outgroups'),
            ('Universal', 'universal')
        ]
        
        for display_name, key in categories:
            stats = summary['length_stats'][key]
            f.write(f"{display_name:22s} {stats['count']:5d}    {stats['mean']:8.0f}   {stats['median']:9.0f}   {stats['min']:7.0f}   {stats['max']:8.0f}   {stats['std']:7.0f}\n")
        
        f.write("\n--- Selection peak length distribution by geographic pattern ---\n")
        f.write("Category                   Count    Mean(bp)   Median(bp)   Min(bp)    Max(bp)    Std(bp)\n")
        f.write("-" * 85 + "\n")
        
        for display_name, key in categories:
            stats = summary['peak_length_stats'][key]
            f.write(f"{display_name:22s} {stats['count']:5d}    {stats['mean']:8.0f}   {stats['median']:9.0f}   {stats['min']:7.0f}   {stats['max']:8.0f}   {stats['std']:7.0f}\n")
        
        f.write("\n--- Present regions per population label ---\n")
        for label, count in sorted(summary['present_per_label'].items()):
            f.write(f"{label:30s}: {count:4d} ({count/summary['total_regions']*100:.1f}%)\n")
        
        f.write("\n--- Weakly present regions per population label ---\n")
        for label, count in sorted(summary['weak_present_per_label'].items()):
            f.write(f"{label:30s}: {count:4d} ({count/summary['total_regions']*100:.1f}%)\n")
        
        f.write("\n--- Average selection regions per population within each label ---\n")
        f.write("Population Label               Populations  Mean±SD      Median   Min   Max\n")
        f.write("-" * 75 + "\n")
        for label in sorted(summary['regions_per_population_stats'].keys()):
            stats_data = summary['regions_per_population_stats'][label]
            f.write(f"{label:30s} {stats_data['num_populations']:11d}  {stats_data['mean_regions_per_pop']:4.1f}±{stats_data['std_regions_per_pop']:4.1f}   {stats_data['median_regions_per_pop']:6.1f}  {stats_data['min_regions_per_pop']:4.0f}  {stats_data['max_regions_per_pop']:4.0f}\n")
        
        # Add Mann-Whitney U test results
        if summary['mann_whitney_test'] is not None:
            test = summary['mann_whitney_test']
            f.write("\n--- Mann-Whitney U Test: North-West vs South-East ---\n")
            f.write("Statistical test for difference in selection regions per population\n")
            f.write(f"North-West populations: n={test['north_west_n']:2d}, mean={test['north_west_mean']:5.1f}, median={test['north_west_median']:5.1f}\n")
            f.write(f"South-East populations: n={test['south_east_n']:2d}, mean={test['south_east_mean']:5.1f}, median={test['south_east_median']:5.1f}\n")
            f.write(f"U-statistic: {test['u_statistic']:8.1f}\n")
            f.write(f"P-value: {test['p_value']:12.6f}")
            if test['p_value'] < 0.001:
                f.write(" (highly significant, p < 0.001)")
            elif test['p_value'] < 0.01:
                f.write(" (highly significant, p < 0.01)")
            elif test['p_value'] < 0.05:
                f.write(" (significant, p < 0.05)")
            else:
                f.write(" (not significant, p ≥ 0.05)")
            f.write("\n")
        else:
            f.write("\n--- Mann-Whitney U Test: North-West vs South-East (All Regions) ---\n")
            f.write("Test not performed: insufficient data in one or both groups\n")
        
        # Add region-specific Mann-Whitney U test results
        if summary['mann_whitney_test_specific'] is not None:
            test = summary['mann_whitney_test_specific']
            f.write("\n--- Mann-Whitney U Test: Region-Specific Selections Only ---\n")
            f.write("Statistical test comparing region-specific selection burden\n")
            f.write(f"Total North-West only regions: {test['north_west_only_regions']}\n")
            f.write(f"Total South-East only regions: {test['south_east_only_regions']}\n")
            f.write(f"North-West populations: n={test['north_west_n']:2d}, mean={test['north_west_mean']:5.1f}, median={test['north_west_median']:5.1f}\n")
            f.write(f"South-East populations: n={test['south_east_n']:2d}, mean={test['south_east_mean']:5.1f}, median={test['south_east_median']:5.1f}\n")
            f.write(f"U-statistic: {test['u_statistic']:8.1f}\n")
            f.write(f"P-value: {test['p_value']:12.6f}")
            if test['p_value'] < 0.001:
                f.write(" (highly significant, p < 0.001)")
            elif test['p_value'] < 0.01:
                f.write(" (highly significant, p < 0.01)")
            elif test['p_value'] < 0.05:
                f.write(" (significant, p < 0.05)")
            else:
                f.write(" (not significant, p ≥ 0.05)")
            f.write("\n")
        else:
            f.write("\n--- Mann-Whitney U Test: Region-Specific Selections Only ---\n")
            f.write("Test not performed: insufficient region-specific data in one or both groups\n")
        
        # Add peak length comparison test results
        if summary['mann_whitney_test_peaks'] is not None:
            test = summary['mann_whitney_test_peaks']
            f.write("\n--- Mann-Whitney U Test: Peak Lengths (Shared vs Unique Regions) ---\n")
            f.write("Statistical test comparing peak lengths between shared and geographically unique regions\n")
            f.write(f"Shared NW-SE regions: {test['shared_regions_count']} regions contributing {test['shared_n']} peaks\n")
            f.write(f"Unique (NW-only + SE-only): {test['unique_regions_count']} regions contributing {test['unique_n']} peaks\n")
            f.write(f"Shared peaks:    mean={test['shared_mean']:6.1f} bp, median={test['shared_median']:6.1f} bp\n")
            f.write(f"Unique peaks:    mean={test['unique_mean']:6.1f} bp, median={test['unique_median']:6.1f} bp\n")
            f.write(f"U-statistic: {test['u_statistic']:8.1f}\n")
            f.write(f"P-value: {test['p_value']:12.6f}")
            if test['p_value'] < 0.001:
                f.write(" (highly significant, p < 0.001)")
            elif test['p_value'] < 0.01:
                f.write(" (highly significant, p < 0.01)")
            elif test['p_value'] < 0.05:
                f.write(" (significant, p < 0.05)")
            else:
                f.write(" (not significant, p ≥ 0.05)")
            f.write("\n")
            
            # Add interpretation hint
            if test['p_value'] < 0.05:
                if test['shared_median'] > test['unique_median']:
                    f.write("Interpretation: Shared regions have significantly longer peaks than unique regions\n")
                else:
                    f.write("Interpretation: Unique regions have significantly longer peaks than shared regions\n")
        else:
            f.write("\n--- Mann-Whitney U Test: Peak Lengths (Shared vs Unique Regions) ---\n")
            f.write("Test not performed: insufficient peak data in one or both groups\n")
        
        f.write("\n--- Detailed distribution patterns ---\n")
        f.write("Format: NW:X_SE:Y_OUT:Z means X North-West labels, Y South-East labels, Z Outgroup labels\n")
        for pattern, count in sorted(summary['detailed_patterns'].items(), key=lambda x: x[1], reverse=True):
            f.write(f"{pattern:20s}: {count:4d} ({count/summary['total_regions']*100:.1f}%)\n")

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Analyze selection regions across population labels')
    parser.add_argument('--threshold', '-t', type=float, default=0.5,
                       help='Overlap fraction threshold for determining presence (default: 0.5)')
    parser.add_argument('--min-length', type=int, default=5000,
                       help='Minimum region length to include in analysis (default: 5000)')
    args = parser.parse_args()
    
    # File paths
    label_file = "data/population_label.txt"
    selection_file = "data/EntNic_grassland_selection_regions.tsv"
    
    # Create output filenames with threshold info
    threshold_str = str(args.threshold).replace('.', 'p')
    output_file_wide = f"data/selection_regions_presence_absence_wide_t{threshold_str}.tsv"
    output_file_long = f"data/selection_regions_presence_absence_long_t{threshold_str}.tsv"
    output_file_summary = f"data/selection_regions_distribution_summary_t{threshold_str}.txt"
    
    # Parse population labels
    population_labels = parse_population_labels(label_file)
    print("Population label mapping:")
    for pop, label in sorted(population_labels.items()):
        print(f"  {pop} -> {label}")
    print()
    
    print(f"Using overlap fraction threshold: {args.threshold}")
    print(f"Using minimum region length: {args.min_length} bp")
    print()
    
    # Analyze selection regions
    results_df, original_filtered_df = analyze_selection_regions(selection_file, population_labels, min_length=args.min_length)
    
    # Categorize presence/absence based on thresholds
    categorized_df = categorize_presence_absence(results_df, threshold_present=args.threshold, threshold_absent=0.1)
    
    # Save wide format results (both fine-grained and categorized)
    results_df.to_csv(output_file_wide, sep='\t', index=False)
    print(f"Wide format results saved to {output_file_wide}")
    
    # Save categorized results
    output_file_categorized = f"data/selection_regions_categorized_t{threshold_str}.tsv"
    categorized_df.to_csv(output_file_categorized, sep='\t', index=False)
    print(f"Categorized results saved to {output_file_categorized}")
    
    # Create and save long format results
    long_df = create_long_format(results_df, population_labels)
    long_df.to_csv(output_file_long, sep='\t', index=False)
    print(f"Long format results saved to {output_file_long}")
    
    # Analyze and write distribution summary
    summary = summarize_selection_distribution(categorized_df, population_labels, original_filtered_df, overlap_threshold=args.threshold, include_weak=True)
    write_distribution_summary(summary, output_file_summary, overlap_threshold=args.threshold)
    print(f"Distribution summary saved to {output_file_summary}")
    
    # Print basic summary
    print(f"\nBasic Summary:")
    print(f"- Total regions >= {args.min_length}bp: {len(results_df)}")
    print(f"- Population labels: {sorted(set(population_labels.values()))}")
    print(f"- Wide format dimensions: {results_df.shape}")
    print(f"- Long format dimensions: {long_df.shape}")
    
    # Show first few rows of categorized results
    print(f"\nFirst 3 regions (categorized):")
    basic_cols = ['region_id', 'chromosome', 'start', 'end', 'length']
    category_cols = [col for col in categorized_df.columns if col.endswith('_category')]
    print(categorized_df[basic_cols + category_cols[:3]].head(3))

if __name__ == "__main__":
    main()
