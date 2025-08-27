#!/usr/bin/env python3
"""
Bootstrap phylogeny sequence sampler for pooled VCF data.
Samples one concatenated sequence per population from all genome-wide SNPs.
"""

import argparse
import random
import sys
import gzip
import os

def parse_vcf_header(vcf_file):
    """Parse VCF header to extract sample names."""
    if vcf_file.endswith('.gz'):
        f = gzip.open(vcf_file, 'rt')
    else:
        f = open(vcf_file, 'r')
    
    samples = []
    for line in f:
        if line.startswith('#CHROM'):
            parts = line.strip().split('\t')
            samples = parts[9:]  # Sample names start from column 9
            break
        if not line.startswith('#'):
            break
    
    f.close()
    return samples

def parse_ad_field(ad_string):
    """Parse AD field to get reference and alternate allele counts."""
    if ad_string == '.' or ad_string == '.,.':
        return 0, 0
    
    parts = ad_string.split(',')
    if len(parts) >= 2:
        try:
            ref_count = int(parts[0]) if parts[0] != '.' else 0
            alt_count = int(parts[1]) if parts[1] != '.' else 0
            return ref_count, alt_count
        except ValueError:
            return 0, 0
    return 0, 0

def sample_allele(ref_count, alt_count, ref_allele, alt_allele):
    """Sample an allele based on read counts (allele frequencies)."""
    total = ref_count + alt_count
    if total == 0:
        return ref_allele  # No coverage means no variants, return reference
    
    # Calculate probability of reference allele
    ref_prob = ref_count / total
    
    # Sample based on probability
    if random.random() < ref_prob:
        return ref_allele
    else:
        return alt_allele

def parse_region(region_str):
    """Parse region string in format 'scaffold:start-end' or just 'scaffold'."""
    if ':' in region_str:
        scaffold, coords = region_str.split(':', 1)
        if '-' in coords:
            start, end = coords.split('-', 1)
            return scaffold, int(start), int(end)
        else:
            return scaffold, int(coords), int(coords)
    else:
        return region_str, None, None

def passes_region_filter(chrom, pos, scaffolds, regions):
    """Check if variant passes scaffold/region filters."""
    # If specific scaffolds specified, check if this chromosome is included
    if scaffolds and chrom not in scaffolds:
        return False
    
    # If specific regions specified, check if position falls within any region
    if regions:
        for region_scaffold, start, end in regions:
            if chrom == region_scaffold:
                if start is None and end is None:
                    return True  # Whole scaffold
                elif start is not None and end is not None:
                    return start <= pos <= end
                elif start is not None:
                    return pos >= start
                elif end is not None:
                    return pos <= end
        return False
    
    return True

def process_vcf(vcf_file, output_file, bootstrap_replicates=1, scaffolds=None, regions=None, outdir=None):
    """Process VCF file and generate bootstrap FASTA file with concatenated sequences."""
    
    print(f"Parsing samples from {vcf_file}...")
    samples = parse_vcf_header(vcf_file)
    print(f"Found {len(samples)} samples")
    
    if scaffolds:
        print(f"Filtering to scaffolds: {scaffolds}")
    if regions:
        print(f"Filtering to regions: {[f'{s}:{start}-{end}' if start and end else s for s, start, end in regions]}")
    
    # Generate bootstrap replicates
    for rep in range(bootstrap_replicates):
        print(f"\nGenerating bootstrap replicate {rep + 1}...")
        
        if vcf_file.endswith('.gz'):
            f = gzip.open(vcf_file, 'rt')
        else:
            f = open(vcf_file, 'r')
        
        # Initialize sequences for each sample
        sequences = {sample: [] for sample in samples}
        
        print("Reading and sampling variants...")
        variant_count = 0
        filtered_count = 0
        
        for line in f:
            if line.startswith('#'):
                continue
                
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue
                
            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            
            # Apply region/scaffold filters
            if not passes_region_filter(chrom, pos, scaffolds, regions):
                filtered_count += 1
                continue
            
            # Skip multi-allelic sites (keep it simple)
            if ',' in alt:
                continue
                
            # Extract genotype data
            format_fields = parts[8].split(':')
            try:
                ad_index = format_fields.index('AD')
            except ValueError:
                continue  # Skip if no AD field
                
            # Process each sample for this variant
            for i, sample in enumerate(samples):
                sample_data = parts[9 + i].split(':')
                if len(sample_data) > ad_index:
                    ref_count, alt_count = parse_ad_field(sample_data[ad_index])
                else:
                    ref_count, alt_count = 0, 0
                
                # Sample allele for this population at this site
                sampled_allele = sample_allele(ref_count, alt_count, ref, alt)
                sequences[sample].append(sampled_allele)
            
            variant_count += 1
            if variant_count % 10000 == 0:
                print(f"  Processed {variant_count} variants (filtered {filtered_count})...")
        
        f.close()
        
        print(f"Total variants processed: {variant_count} (filtered out: {filtered_count})")
        
        # Write output FASTA file
        if bootstrap_replicates > 1:
            out_filename = f"{output_file.rsplit('.', 1)[0]}_rep{rep+1:03d}.fasta"
        else:
            out_filename = output_file if output_file.endswith('.fasta') else f"{output_file}.fasta"
        
        # Add output directory if specified
        if outdir:
            out_filename = os.path.join(outdir, os.path.basename(out_filename))
        
        print(f"Writing sequences to {out_filename}...")
        
        with open(out_filename, 'w') as out_f:
            for sample in samples:
                # Write FASTA header
                header = f">{sample}"
                if bootstrap_replicates > 1:
                    header += f"_rep{rep+1:03d}"
                out_f.write(header + '\n')
                
                # Write concatenated sequence in 80-character lines
                seq_str = ''.join(sequences[sample])
                for i in range(0, len(seq_str), 80):
                    out_f.write(seq_str[i:i+80] + '\n')
        
        print(f"Written {len(samples)} concatenated sequences of length {len(sequences[samples[0]])} to {out_filename}")

def main():
    parser = argparse.ArgumentParser(description='Sample bootstrap sequences for phylogeny from pooled VCF')
    parser.add_argument('vcf', help='Input VCF file (can be gzipped)')
    parser.add_argument('-o', '--output', default='phylo_sequences.fasta', 
                       help='Output FASTA file (default: phylo_sequences.fasta)')
    parser.add_argument('-r', '--replicates', type=int, default=1,
                       help='Number of bootstrap replicates (default: 1)')
    parser.add_argument('-s', '--scaffolds', nargs='+',
                       help='Specific scaffolds to include (e.g., HiC_scaffold_1 HiC_scaffold_2)')
    parser.add_argument('--regions', nargs='+',
                       help='Specific regions to include (format: scaffold:start-end or scaffold)')
    parser.add_argument('--seed', type=int, help='Random seed for reproducibility')
    parser.add_argument('--outdir', help='Output directory for FASTA files')
    
    args = parser.parse_args()
    
    if args.seed:
        random.seed(args.seed)
        print(f"Random seed set to: {args.seed}")
    
    # Create output directory if specified
    if args.outdir:
        os.makedirs(args.outdir, exist_ok=True)
        print(f"Output directory: {args.outdir}")
    
    # Parse regions if provided
    regions = None
    if args.regions:
        regions = []
        for region_str in args.regions:
            regions.append(parse_region(region_str))
    
    process_vcf(args.vcf, args.output, args.replicates, args.scaffolds, regions, args.outdir)
    print("\nBootstrap sampling completed!")

if __name__ == '__main__':
    main()