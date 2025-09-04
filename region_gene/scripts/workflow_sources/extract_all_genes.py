#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging
import math

def setup_logging(log_file=None):
    """Set up logging configuration"""
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    if log_file:
        logging.basicConfig(level=logging.INFO, format=log_format, filename=log_file, filemode='w')
    else:
        logging.basicConfig(level=logging.INFO, format=log_format)

def parse_gff_all_genes(gff_file):
    """Parse GFF file to extract all genes in the genome"""
    genes = []
    
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
                
                # Only process genes
                if feature_type.lower() != 'gene':
                    continue
                
                # Extract gene ID from attributes
                gene_id = None
                for attr in attributes.split(';'):
                    if attr.startswith('ID='):
                        gene_id = attr.split('=')[1]
                        break
                
                if gene_id:
                    genes.append({
                        'gene_id': gene_id,
                        'start': feature_start,
                        'end': feature_end,
                        'strand': strand,
                        'chromosome': seqid
                    })
    
    except Exception as e:
        logging.error(f"Error parsing GFF file {gff_file}: {e}")
        return []
    
    logging.info(f"Found {len(genes)} genes in total across all chromosomes")
    return genes

def extract_gene_sequences(genome_file, genes, output_file):
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
            
            # Create systematic sequence record name
            seq_record = SeqRecord(
                gene_seq,
                id=f"{gene_id}_{chromosome}_{gene['start']}_{gene['end']}",
                description=f"Gene:{gene_id} Coords:{chromosome}:{gene['start']}-{gene['end']} Strand:{strand} Length:{len(gene_seq)}"
            )
            
            sequences.append(seq_record)
        
        # Write sequences to file
        if sequences:
            SeqIO.write(sequences, output_file, "fasta")
            logging.info(f"Wrote {len(sequences)} gene sequences to {output_file}")
            return len(sequences)
        else:
            logging.warning("No gene sequences found")
            # Create empty file to indicate completion
            with open(output_file, 'w') as f:
                f.write("# No genes found in genome\n")
            return 0
    
    except Exception as e:
        logging.error(f"Error extracting gene sequences: {e}")
        return -1

def create_gene_chunks(genes, chunk_size):
    """Split genes into chunks for parallel processing"""
    chunks = []
    for i in range(0, len(genes), chunk_size):
        chunk = genes[i:i + chunk_size]
        chunk_id = f"{i//chunk_size + 1:04d}"
        chunks.append((chunk_id, chunk))
    
    logging.info(f"Created {len(chunks)} chunks with chunk size {chunk_size}")
    return chunks

def main():
    parser = argparse.ArgumentParser(description='Extract all gene sequences from genome using GFF annotation')
    parser.add_argument('--genome-file', required=True, help='Genome FASTA file')
    parser.add_argument('--gff-file', required=True, help='Genome annotation GFF file')
    parser.add_argument('--output-dir', required=True, help='Output directory for FASTA files')
    parser.add_argument('--log-file', help='Log file path')
    parser.add_argument('--output-filename', default='all_genes.fasta', help='Output filename for all genes (default: all_genes.fasta)')
    
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(args.log_file)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    try:
        # Read all genes from GFF file
        logging.info(f"Extracting all genes from {args.gff_file}")
        all_genes = parse_gff_all_genes(args.gff_file)
        
        if not all_genes:
            logging.error("No genes found in GFF file")
            sys.exit(1)
        
        # Extract all gene sequences into a single file
        output_file = os.path.join(args.output_dir, args.output_filename)
        logging.info(f"Processing all {len(all_genes)} genes into single file: {output_file}")
        
        num_genes = extract_gene_sequences(args.genome_file, all_genes, output_file)
        
        if num_genes >= 0:
            # Create summary file
            summary_file = os.path.join(args.output_dir, "extraction_summary.txt")
            with open(summary_file, 'w') as f:
                f.write(f"Total genes in GFF: {len(all_genes)}\n")
                f.write(f"Total unique genes processed: {num_genes}\n")
                f.write(f"Output file: {output_file}\n")
            
            logging.info(f"Successfully processed {num_genes} unique genes from {len(all_genes)} total genes")
            logging.info(f"Output written to: {output_file}")
            logging.info(f"Summary written to: {summary_file}")
        else:
            logging.error("Failed to process genes")
            sys.exit(1)
    
    except Exception as e:
        logging.error(f"Error in main processing: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()