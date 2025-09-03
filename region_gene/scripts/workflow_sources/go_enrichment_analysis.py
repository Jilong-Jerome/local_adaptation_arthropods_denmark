#!/usr/bin/env python3

import argparse
import pandas as pd
import requests
import time
import sys
from collections import defaultdict
import logging

def setup_logging(log_file=None):
    """Set up logging configuration"""
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    if log_file:
        logging.basicConfig(level=logging.INFO, format=log_format, filename=log_file, filemode='w')
    else:
        logging.basicConfig(level=logging.INFO, format=log_format)

def get_go_terms_batch(uniprot_ids, batch_size=100):
    """Retrieve GO terms for UniProt IDs in batches"""
    go_data = []
    
    for i in range(0, len(uniprot_ids), batch_size):
        batch = uniprot_ids[i:i + batch_size]
        logging.info(f"Processing batch {i//batch_size + 1}: {len(batch)} IDs")
        
        # Create query string
        query = ' OR '.join([f'accession:{uid}' for uid in batch])
        
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            'query': query,
            'format': 'json',
            'fields': 'accession,go_p,go_c,go_f'  # GO Process, Component, Function
        }
        
        try:
            response = requests.get(url, params=params, timeout=30)
            if response.status_code == 200:
                data = response.json()
                
                for entry in data.get('results', []):
                    accession = entry.get('primaryAccession', '')
                    
                    # Process GO terms
                    for go_type in ['goAnnotations']:
                        if go_type in entry:
                            for go_term in entry[go_type]:
                                go_data.append({
                                    'uniprot_id': accession,
                                    'go_id': go_term.get('id', ''),
                                    'go_term': go_term.get('properties', [{}])[0].get('term', '') if go_term.get('properties') else '',
                                    'go_aspect': go_term.get('aspect', ''),
                                    'evidence': go_term.get('evidenceCode', '')
                                })
                
                logging.info(f"Retrieved GO terms for batch {i//batch_size + 1}")
                
            else:
                logging.warning(f"Failed to retrieve batch {i//batch_size + 1}: {response.status_code}")
        
        except Exception as e:
            logging.error(f"Error processing batch {i//batch_size + 1}: {e}")
        
        # Rate limiting
        time.sleep(1)
    
    return pd.DataFrame(go_data)

def perform_go_enrichment(go_df, background_size=20000):
    """Simple GO term enrichment analysis"""
    
    # Count GO terms in our dataset
    go_counts = go_df.groupby(['go_id', 'go_term', 'go_aspect']).size().reset_index(name='count')
    
    # Calculate basic enrichment metrics
    total_proteins = go_df['uniprot_id'].nunique()
    
    enrichment_results = []
    
    for _, row in go_counts.iterrows():
        go_id = row['go_id']
        go_term = row['go_term']
        go_aspect = row['go_aspect']
        count = row['count']
        
        # Simple enrichment calculation
        frequency = count / total_proteins
        
        # Basic statistical significance (simplified)
        # In a real analysis, you'd use hypergeometric test or Fisher's exact test
        enrichment_score = frequency * 100  # Convert to percentage
        
        enrichment_results.append({
            'GO_ID': go_id,
            'GO_Term': go_term,
            'GO_Aspect': go_aspect,
            'Count': count,
            'Frequency': frequency,
            'Percentage': enrichment_score,
            'Total_Proteins': total_proteins
        })
    
    return pd.DataFrame(enrichment_results).sort_values('Count', ascending=False)

def main():
    parser = argparse.ArgumentParser(description='GO term enrichment analysis for UniProt IDs')
    parser.add_argument('--input-file', required=True, help='Input file with UniProt IDs (one per line)')
    parser.add_argument('--output-dir', required=True, help='Output directory for results')
    parser.add_argument('--log-file', help='Log file path')
    parser.add_argument('--batch-size', type=int, default=50, help='Batch size for API requests (default: 50)')
    
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(args.log_file)
    
    # Create output directory
    import os
    os.makedirs(args.output_dir, exist_ok=True)
    
    try:
        # Read UniProt IDs
        logging.info(f"Reading UniProt IDs from {args.input_file}")
        with open(args.input_file, 'r') as f:
            uniprot_ids = [line.strip() for line in f if line.strip()]
        
        # Remove duplicates while preserving order
        unique_ids = list(dict.fromkeys(uniprot_ids))
        
        logging.info(f"Found {len(uniprot_ids)} total IDs, {len(unique_ids)} unique IDs")
        
        # Retrieve GO terms
        logging.info("Starting GO term retrieval...")
        go_df = get_go_terms_batch(unique_ids, args.batch_size)
        
        if go_df.empty:
            logging.warning("No GO terms retrieved. Check your UniProt IDs and internet connection.")
            return
        
        logging.info(f"Retrieved {len(go_df)} GO term associations")
        
        # Save raw GO data
        go_output = os.path.join(args.output_dir, 'go_annotations.tsv')
        go_df.to_csv(go_output, sep='\t', index=False)
        logging.info(f"Raw GO annotations saved to {go_output}")
        
        # Perform enrichment analysis
        logging.info("Performing GO term enrichment analysis...")
        enrichment_df = perform_go_enrichment(go_df)
        
        # Save enrichment results
        enrichment_output = os.path.join(args.output_dir, 'go_enrichment_results.tsv')
        enrichment_df.to_csv(enrichment_output, sep='\t', index=False)
        logging.info(f"GO enrichment results saved to {enrichment_output}")
        
        # Summary by GO aspect
        summary_by_aspect = go_df.groupby('go_aspect').agg({
            'uniprot_id': 'nunique',
            'go_id': 'nunique'
        }).rename(columns={'uniprot_id': 'Unique_Proteins', 'go_id': 'Unique_GO_Terms'})
        
        summary_output = os.path.join(args.output_dir, 'go_summary_by_aspect.tsv')
        summary_by_aspect.to_csv(summary_output, sep='\t')
        logging.info(f"GO summary by aspect saved to {summary_output}")
        
        # Top enriched terms by aspect
        for aspect in ['P', 'F', 'C']:  # Process, Function, Component
            aspect_data = enrichment_df[enrichment_df['GO_Aspect'] == aspect].head(20)
            if not aspect_data.empty:
                aspect_names = {'P': 'Biological_Process', 'F': 'Molecular_Function', 'C': 'Cellular_Component'}
                aspect_output = os.path.join(args.output_dir, f'top_20_{aspect_names[aspect]}_terms.tsv')
                aspect_data.to_csv(aspect_output, sep='\t', index=False)
                logging.info(f"Top 20 {aspect_names[aspect]} terms saved to {aspect_output}")
        
        # Final summary
        total_proteins_with_go = go_df['uniprot_id'].nunique()
        total_go_terms = go_df['go_id'].nunique()
        total_associations = len(go_df)
        
        logging.info("="*50)
        logging.info("GO ENRICHMENT ANALYSIS SUMMARY")
        logging.info("="*50)
        logging.info(f"Input UniProt IDs: {len(unique_ids)}")
        logging.info(f"Proteins with GO annotations: {total_proteins_with_go}")
        logging.info(f"Unique GO terms found: {total_go_terms}")
        logging.info(f"Total GO associations: {total_associations}")
        logging.info(f"Coverage: {total_proteins_with_go/len(unique_ids)*100:.1f}%")
        
    except Exception as e:
        logging.error(f"Error in GO enrichment analysis: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()