#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging

def setup_logging(log_file=None):
    """Set up logging configuration"""
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    if log_file:
        logging.basicConfig(level=logging.INFO, format=log_format, filename=log_file, filemode='w')
    else:
        logging.basicConfig(level=logging.INFO, format=log_format)

def parse_gff_genes(gff_file, chromosome, start, end, flank_size=0):
    """Parse GFF file to find genes within a specific region (with optional flanking)"""
    genes = []
    
    # Extend search region with flanking sequences
    extended_start = max(1, start - flank_size)  # Don't go below 1
    extended_end = end + flank_size
    
    try:
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                seqid = parts[0]
                feature_type = parts[2]
                feature_start = int(parts[3])
                feature_end = int(parts[4])
                strand = parts[6]
                attributes = parts[8]
                
                # Only process genes on the target chromosome
                if seqid != chromosome or feature_type.lower() != 'gene':
                    continue
                
                # Check if gene overlaps with extended region
                if feature_start <= extended_end and feature_end >= extended_start:
                    # Extract gene ID from attributes
                    gene_id = None
                    for attr in attributes.split(';'):
                        if attr.startswith('ID='):
                            gene_id = attr.split('=')[1]
                            break
                    
                    if gene_id:
                        # Determine gene position relative to original region
                        if feature_end < start:
                            position = "upstream"
                        elif feature_start > end:
                            position = "downstream"  
                        else:
                            position = "overlapping"
                        
                        genes.append({
                            'gene_id': gene_id,
                            'start': feature_start,
                            'end': feature_end,
                            'strand': strand,
                            'chromosome': seqid,
                            'position': position
                        })
    
    except Exception as e:
        logging.error(f"Error parsing GFF file {gff_file}: {e}")
        return []
    
    return genes

def extract_gene_sequences(genome_file, genes, output_file, region_id, region_start, region_end):
    """Extract gene sequences from genome and write to FASTA file"""
    try:
        # Load genome sequences
        genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
        logging.info(f"Loaded {len(genome)} sequences from genome file")
        
        sequences = []
        seen_locations = set()  # Track duplicates by genomic location only
        
        for gene in genes:
            chromosome = gene['chromosome']
            gene_id = gene['gene_id']
            position = gene['position']
            
            # Check for duplicates (same genomic coordinates only)
            location_key = f"{chromosome}_{gene['start']}_{gene['end']}"
            if location_key in seen_locations:
                logging.warning(f"Duplicate location found: {chromosome}:{gene['start']}-{gene['end']} (gene_id: {gene_id})")
                continue
            seen_locations.add(location_key)
            
            if chromosome not in genome:
                logging.warning(f"Chromosome {chromosome} not found in genome file")
                continue
            
            chrom_seq = genome[chromosome]
            gene_start = gene['start'] - 1  # Convert to 0-based indexing
            gene_end = gene['end']
            strand = gene['strand']
            
            # Extract sequence
            gene_seq = chrom_seq.seq[gene_start:gene_end]
            
            # Reverse complement if on negative strand
            if strand == '-':
                gene_seq = gene_seq.reverse_complement()
            
            # Calculate region length
            region_length = region_end - region_start + 1
            
            # Create systematic sequence record name
            seq_record = SeqRecord(
                gene_seq,
                id=f"{region_id}_R{region_start}-{region_end}_L{region_length}_{position}_{gene_id}_{chromosome}_{gene['start']}_{gene['end']}",
                description=f"Region:{region_id}({region_start}-{region_end},L{region_length}) Position:{position} Gene:{gene_id} Coords:{chromosome}:{gene['start']}-{gene['end']} Strand:{strand}"
            )
            
            sequences.append(seq_record)
        
        # Write sequences to file
        if sequences:
            SeqIO.write(sequences, output_file, "fasta")
            logging.info(f"Wrote {len(sequences)} gene sequences to {output_file}")
            return len(sequences)
        else:
            logging.warning(f"No gene sequences found for region {region_id}")
            # Create empty file to indicate completion
            with open(output_file, 'w') as f:
                f.write(f"# No genes found in region {region_id}\n")
            return 0
    
    except Exception as e:
        logging.error(f"Error extracting gene sequences: {e}")
        return -1

def main():
    parser = argparse.ArgumentParser(description='Extract gene sequences from selection regions')
    parser.add_argument('--regions-file', required=True, help='TSV file with selection regions')
    parser.add_argument('--genome-file', required=True, help='Genome FASTA file')
    parser.add_argument('--gff-file', required=True, help='Genome annotation GFF file')
    parser.add_argument('--output-dir', required=True, help='Output directory for FASTA files')
    parser.add_argument('--region-id', help='Specific region ID to process (if not provided, process all)')
    parser.add_argument('--log-file', help='Log file path')
    parser.add_argument('--min-length', type=int, default=5000, help='Minimum region length threshold (default: 5000)')
    parser.add_argument('--flank-size', type=int, default=10000, help='Flanking region size in bp (default: 10000)')
    
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(args.log_file)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    try:
        # Read regions file
        logging.info(f"Reading regions from {args.regions_file}")
        regions_df = pd.read_csv(args.regions_file, sep='\t')
        
        # Get unique regions
        unique_regions = regions_df[['region_id', 'chromosome', 'start', 'end']].drop_duplicates()
        
        # Filter by minimum length
        unique_regions['length'] = unique_regions['end'] - unique_regions['start'] + 1
        before_filter = len(unique_regions)
        unique_regions = unique_regions[unique_regions['length'] >= args.min_length]
        after_filter = len(unique_regions)
        logging.info(f"Filtered regions by minimum length {args.min_length}bp: {before_filter} -> {after_filter}")
        
        if args.region_id:
            # Process specific region
            region_row = unique_regions[unique_regions['region_id'] == args.region_id]
            if region_row.empty:
                logging.error(f"Region {args.region_id} not found in regions file or below minimum length threshold")
                sys.exit(1)
            regions_to_process = region_row
        else:
            # Process all regions
            regions_to_process = unique_regions
        
        total_processed = 0
        total_genes = 0
        
        for _, region in regions_to_process.iterrows():
            region_id = region['region_id']
            chromosome = region['chromosome']
            start = int(region['start'])
            end = int(region['end'])
            
            logging.info(f"Processing region {region_id}: {chromosome}:{start}-{end} (length: {end-start+1}bp, +/-{args.flank_size}bp flanks)")
            
            # Find genes in this region (with flanking)
            genes = parse_gff_genes(args.gff_file, chromosome, start, end, args.flank_size)
            logging.info(f"Found {len(genes)} genes in region {region_id}")
            
            # Extract gene sequences
            output_file = os.path.join(args.output_dir, f"{region_id}_genes.fasta")
            num_genes = extract_gene_sequences(args.genome_file, genes, output_file, region_id, start, end)
            
            if num_genes >= 0:
                total_processed += 1
                total_genes += num_genes
            else:
                logging.error(f"Failed to process region {region_id}")
        
        logging.info(f"Successfully processed {total_processed} regions with a total of {total_genes} genes")
    
    except Exception as e:
        logging.error(f"Error in main processing: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()