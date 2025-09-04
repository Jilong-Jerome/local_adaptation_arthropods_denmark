#!/usr/bin/env python3
"""
Calculate total size of selection regions with optional flanking regions.
"""

import pandas as pd
import argparse


def calculate_region_sizes(tsv_file, flanking_size=0, min_size=0):
    """
    Calculate total size of unique selection regions.
    
    Args:
        tsv_file: Path to the TSV file containing selection regions
        flanking_size: Size of flanking regions to add on both ends (default: 0)
        min_size: Minimum size threshold for selection regions (default: 0)
    
    Returns:
        tuple: (total_original_size, total_with_flanking_size, num_regions)
    """
    # Read the TSV file
    df = pd.read_csv(tsv_file, sep='\t')
    
    # Get unique regions based on region_id, chromosome, start, end
    unique_regions = df[['region_id', 'chromosome', 'start', 'end']].drop_duplicates()
    
    # Calculate original region sizes
    unique_regions['size'] = unique_regions['end'] - unique_regions['start'] + 1
    
    # Filter by minimum size threshold
    filtered_regions = unique_regions[unique_regions['size'] >= min_size].copy()
    total_original_size = filtered_regions['size'].sum()
    
    # Calculate sizes with flanking regions
    filtered_regions['size_with_flanking'] = (
        filtered_regions['end'] + flanking_size - 
        (filtered_regions['start'] - flanking_size) + 1
    )
    total_with_flanking_size = filtered_regions['size_with_flanking'].sum()
    
    return total_original_size, total_with_flanking_size, len(filtered_regions), len(unique_regions)


def main():
    parser = argparse.ArgumentParser(
        description="Calculate total size of selection regions with optional flanking regions"
    )
    parser.add_argument(
        "tsv_file", 
        help="Path to the TSV file containing selection regions"
    )
    parser.add_argument(
        "--flanking", 
        type=int, 
        default=10000,
        help="Size of flanking regions to add on both ends (default: 10000)"
    )
    parser.add_argument(
        "--min-size", 
        type=int, 
        default=5000,
        help="Minimum size threshold for selection regions in bp (default: 5000)"
    )
    
    args = parser.parse_args()
    
    try:
        original_size, flanking_size, num_filtered_regions, num_total_regions = calculate_region_sizes(
            args.tsv_file, args.flanking, args.min_size
        )
        
        print(f"Selection Regions Analysis")
        print(f"=" * 50)
        print(f"Total unique selection regions: {num_total_regions:,}")
        print(f"Regions >= {args.min_size:,} bp: {num_filtered_regions:,}")
        print(f"Regions filtered out: {num_total_regions - num_filtered_regions:,}")
        print(f"Total size of filtered regions: {original_size:,} bp")
        print(f"Flanking region size (each end): {args.flanking:,} bp")
        print(f"Total size with flanking regions: {flanking_size:,} bp")
        print(f"Additional size from flanking: {flanking_size - original_size:,} bp")
        if original_size > 0:
            print(f"Fold increase: {flanking_size / original_size:.2f}x")
        else:
            print(f"Fold increase: N/A (no regions passed size filter)")
        
    except Exception as e:
        print(f"Error processing file: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())