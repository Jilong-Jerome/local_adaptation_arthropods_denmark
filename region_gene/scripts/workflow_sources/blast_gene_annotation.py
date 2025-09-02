#!/usr/bin/env python3

import argparse
import pandas as pd
import subprocess
import os
import sys
import logging
from Bio import SeqIO
import tempfile

def setup_logging(log_file=None):
    """Set up logging configuration"""
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    if log_file:
        logging.basicConfig(level=logging.INFO, format=log_format, filename=log_file, filemode='w')
    else:
        logging.basicConfig(level=logging.INFO, format=log_format)

def run_diamond_blastx(fasta_file, diamond_db, output_file, evalue=1e-5, max_target_seqs=1):
    """Run diamond blastx search against UniProt database"""
    try:
        # diamond blastx command
        cmd = [
            'diamond', 'blastx',
            '--query', fasta_file,
            '--db', diamond_db,
            '--out', output_file,
            '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle',
            '--evalue', str(evalue),
            '--max-target-seqs', str(max_target_seqs),
            '--threads', '4'
        ]
        
        logging.info(f"Running Diamond command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logging.error(f"Diamond failed with error: {result.stderr}")
            return False
        
        logging.info(f"Diamond completed successfully, output written to {output_file}")
        return True
        
    except Exception as e:
        logging.error(f"Error running Diamond: {e}")
        return False

def parse_diamond_results(diamond_output_file, fasta_file):
    """Parse Diamond results and create annotation summary"""
    
    # Read input sequences to get all query IDs
    query_sequences = {}
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            query_sequences[record.id] = len(record.seq)
    except Exception as e:
        logging.error(f"Error reading input FASTA: {e}")
        return None
    
    # Initialize results dictionary
    results = []
    
    # Read Diamond results if file exists and has content
    diamond_hits = {}
    if os.path.exists(diamond_output_file) and os.path.getsize(diamond_output_file) > 0:
        try:
            diamond_df = pd.read_csv(diamond_output_file, sep='\t', header=None,
                                   names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                                         'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                                         'evalue', 'bitscore', 'stitle'])
            
            # Group by query sequence and take best hit (lowest evalue, highest bitscore)
            for qseqid in diamond_df['qseqid'].unique():
                query_hits = diamond_df[diamond_df['qseqid'] == qseqid]
                best_hit = query_hits.loc[query_hits['evalue'].idxmin()]
                diamond_hits[qseqid] = best_hit
                
        except Exception as e:
            logging.error(f"Error parsing Diamond results: {e}")
    
    # Process all query sequences
    for query_id, seq_length in query_sequences.items():
        if query_id in diamond_hits:
            hit = diamond_hits[query_id]
            
            # Parse subject title to extract protein information
            stitle = hit['stitle']
            
            # Extract protein information from UniProt format
            # UniProt format: tr|A0A123|GENE_ORGANISM Gene name OS=Organism name OX=taxid GN=gene PE=level SV=version
            organism = "Unknown"
            gene_name = "Unknown"
            protein_name = "Unknown"
            gene_function = stitle
            
            # Parse UniProt format
            if 'OS=' in stitle:
                # Extract organism
                os_start = stitle.find('OS=')
                ox_start = stitle.find('OX=', os_start)
                if os_start != -1 and ox_start != -1:
                    organism = stitle[os_start+3:ox_start].strip()
            
            if 'GN=' in stitle:
                # Extract gene name
                gn_start = stitle.find('GN=')
                pe_start = stitle.find('PE=', gn_start)
                if gn_start != -1 and pe_start != -1:
                    gene_name = stitle[gn_start+3:pe_start].strip()
            
            # Extract protein name (everything before OS=)
            if 'OS=' in stitle:
                protein_name = stitle[:stitle.find('OS=')].strip()
                # Remove UniProt ID prefix if present
                if '|' in protein_name:
                    parts = protein_name.split('|')
                    if len(parts) >= 3:
                        protein_name = parts[2].strip()
            
            gene_function = protein_name
            
            results.append({
                'query_id': query_id,
                'query_length': seq_length,
                'hit_status': 'hit',
                'subject_id': hit['sseqid'],
                'gene_name': gene_name,
                'organism': organism,
                'gene_function': gene_function,
                'percent_identity': hit['pident'],
                'alignment_length': hit['length'],
                'evalue': hit['evalue'],
                'bitscore': hit['bitscore'],
                'query_coverage': round((hit['qend'] - hit['qstart'] + 1) / seq_length * 100, 2)
            })
        else:
            # No hit found
            results.append({
                'query_id': query_id,
                'query_length': seq_length,
                'hit_status': 'no_hit',
                'subject_id': 'N/A',
                'gene_name': 'N/A',
                'organism': 'N/A',
                'gene_function': 'N/A',
                'percent_identity': 'N/A',
                'alignment_length': 'N/A',
                'evalue': 'N/A',
                'bitscore': 'N/A',
                'query_coverage': 'N/A'
            })
    
    return pd.DataFrame(results)

def main():
    parser = argparse.ArgumentParser(description='Diamond search of gene sequences against UniProt database for annotation')
    parser.add_argument('--fasta-file', required=True, help='Input FASTA file with gene sequences')
    parser.add_argument('--diamond-db', required=True, help='Path to Diamond database')
    parser.add_argument('--output-file', required=True, help='Output TSV file with annotation results')
    parser.add_argument('--log-file', help='Log file path')
    parser.add_argument('--evalue', type=float, default=1e-5, help='E-value threshold for Diamond (default: 1e-5)')
    parser.add_argument('--max-hits', type=int, default=1, help='Maximum number of hits per query (default: 1)')
    
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(args.log_file)
    
    # Check if input files exist
    if not os.path.exists(args.fasta_file):
        logging.error(f"Input FASTA file not found: {args.fasta_file}")
        sys.exit(1)
    
    # Check if FASTA file has sequences
    try:
        seq_count = sum(1 for _ in SeqIO.parse(args.fasta_file, "fasta"))
        if seq_count == 0:
            logging.warning(f"No sequences found in {args.fasta_file}")
            # Create empty results file
            empty_df = pd.DataFrame(columns=['query_id', 'query_length', 'hit_status', 'subject_id', 
                                           'gene_name', 'organism', 'gene_function', 'percent_identity', 
                                           'alignment_length', 'evalue', 'bitscore', 'query_coverage'])
            empty_df.to_csv(args.output_file, sep='\t', index=False)
            return
        
        logging.info(f"Found {seq_count} sequences in input FASTA")
        
    except Exception as e:
        logging.error(f"Error reading input FASTA: {e}")
        sys.exit(1)
    
    # Create temporary file for BLAST output
    with tempfile.NamedTemporaryFile(mode='w', suffix='.blast', delete=False) as tmp_file:
        blast_output = tmp_file.name
    
    try:
        # Run Diamond search
        logging.info(f"Starting Diamond search against {args.diamond_db}")
        success = run_diamond_blastx(args.fasta_file, args.diamond_db, blast_output, 
                                   args.evalue, args.max_hits)
        
        if not success:
            logging.error("Diamond search failed")
            sys.exit(1)
        
        # Parse results and create annotation summary
        logging.info("Parsing Diamond results and creating annotation summary")
        results_df = parse_diamond_results(blast_output, args.fasta_file)
        
        if results_df is None:
            logging.error("Failed to parse Diamond results")
            sys.exit(1)
        
        # Write results to output file
        results_df.to_csv(args.output_file, sep='\t', index=False)
        
        # Log summary
        total_queries = len(results_df)
        hits = len(results_df[results_df['hit_status'] == 'hit'])
        no_hits = len(results_df[results_df['hit_status'] == 'no_hit'])
        
        logging.info(f"Diamond annotation completed:")
        logging.info(f"  Total queries: {total_queries}")
        logging.info(f"  Hits found: {hits}")
        logging.info(f"  No hits: {no_hits}")
        logging.info(f"  Results written to: {args.output_file}")
        
    finally:
        # Clean up temporary file
        if os.path.exists(blast_output):
            os.unlink(blast_output)

if __name__ == '__main__':
    main()
