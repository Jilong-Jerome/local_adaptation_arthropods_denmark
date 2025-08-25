#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
from pathlib import Path

def get_population_list(input_file):
    """Get list of unique populations from the file."""
    print("Getting population list...")
    populations = set()
    
    with open(input_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            if line.strip():
                pop = line.split('\t')[4]  # pop column is 5th (index 4)
                populations.add(pop)
    
    populations.discard('pop')  # Remove header if present
    return sorted(list(populations))

def process_population(input_file, population, output_dir):
    """Process likelihood ratio data for a single population."""
    print(f"Processing population: {population}")
    
    lr_values = []
    
    with open(input_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            if line.strip():
                parts = line.split('\t')
                if len(parts) >= 5 and parts[4] == population:
                    lr_values.append(float(parts[1]))  # LR column is 2nd (index 1)
    
    lr_values = np.array(lr_values)
    
    if len(lr_values) == 0:
        print(f"Warning: No data found for population {population}")
        return None
    
    # Generate percentiles from 0.1% to 99.9% with 0.1% resolution
    percentiles = np.arange(0.1, 100.0, 0.1)
    
    # Basic statistics
    stats = {
        'population': population,
        'total_sites': len(lr_values),
        'mean_LR': np.mean(lr_values),
        'median_LR': np.median(lr_values),
        'std_LR': np.std(lr_values),
        'min_LR': np.min(lr_values),
        'max_LR': np.max(lr_values),
        'sites_LR_gt_1': np.sum(lr_values > 1),
        'sites_LR_gt_5': np.sum(lr_values > 5),
        'sites_LR_gt_10': np.sum(lr_values > 10),
        'pct_LR_gt_1': (np.sum(lr_values > 1) / len(lr_values)) * 100,
        'pct_LR_gt_5': (np.sum(lr_values > 5) / len(lr_values)) * 100,
        'pct_LR_gt_10': (np.sum(lr_values > 10) / len(lr_values)) * 100
    }
    
    # Add all percentiles with 0.1% resolution
    percentile_values = np.percentile(lr_values, percentiles)
    for i, perc in enumerate(percentiles):
        stats[f'q{perc:.1f}_LR'] = percentile_values[i]
    
    # Save individual population distribution
    output_dir = Path(output_dir)
    pop_file = output_dir / f"{population}_LR_distribution.txt"
    np.savetxt(pop_file, lr_values, fmt='%.6f', 
              header=f"Likelihood Ratio values for {population}\nTotal sites: {len(lr_values)}")
    
    return stats

def summarize_lr_distribution(input_file, output_dir):
    """
    Summarize likelihood ratio distribution for each population from SweepFinder2 output.
    Process one population at a time to handle large files efficiently.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    print(f"Processing {input_file}...")
    
    # Get list of all populations
    populations = get_population_list(input_file)
    print(f"Found {len(populations)} populations: {populations}")
    
    # Process each population separately
    summary_data = []
    for pop in populations:
        stats = process_population(input_file, pop, output_dir)
        if stats:
            summary_data.append(stats)
    
    # Save summary statistics
    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values('population')
    
    summary_file = output_dir / "LR_summary_statistics.tsv"
    summary_df.to_csv(summary_file, sep='\t', index=False, float_format='%.6f')
    
    print(f"\nSummary completed!")
    print(f"- Summary statistics saved to: {summary_file}")
    print(f"- Individual population distributions saved to: {output_dir}")
    print(f"- Total populations processed: {len(summary_data)}")
    
    return summary_df

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python summarize_lr_distribution.py <input_file> <output_directory>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    
    summary_df = summarize_lr_distribution(input_file, output_dir)
    
    # Print quick overview
    print("\nQuick overview:")
    print(f"Population count: {len(summary_df)}")
    print(f"Mean LR range: {summary_df['mean_LR'].min():.3f} - {summary_df['mean_LR'].max():.3f}")
    print(f"Max LR overall: {summary_df['max_LR'].max():.3f}")
