#!/usr/bin/env python3
"""
GO Enrichment Analysis Script
Performs GO enrichment analysis separating by GO types (C, F, P)
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, false_discovery_control
from collections import defaultdict
import argparse
import sys

def parse_go_terms(go_ids_str, go_terms_str):
    """Parse GO IDs and terms from the TSV format"""
    if pd.isna(go_ids_str) or pd.isna(go_terms_str) or go_ids_str == 'NA':
        return {}
    
    go_dict = {}
    go_ids = go_ids_str.split(';')
    go_terms = go_terms_str.split(';')
    
    for go_id, go_term in zip(go_ids, go_terms):
        # Extract GO type from term (C:, F:, P:)
        if ':' in go_term:
            go_type = go_term[0]  # C, F, or P
            term_name = go_term[2:]  # Remove "C:", "F:", "P:"
            go_dict[go_id.strip()] = {'type': go_type, 'name': term_name}
    
    return go_dict

def load_go_data(filename):
    """Load GO data from TSV file"""
    df = pd.read_csv(filename, sep='\t')
    go_data = {}
    
    for _, row in df.iterrows():
        gene_id = row['Entry']
        if gene_id == 'NA' or pd.isna(gene_id):
            continue
            
        go_terms = parse_go_terms(row['GO_IDs'], row['GO_Terms'])
        if go_terms:
            go_data[gene_id] = go_terms
    
    return go_data

def separate_go_by_type(go_data):
    """Separate GO terms by type (C, F, P)"""
    go_by_type = {'C': defaultdict(set), 'F': defaultdict(set), 'P': defaultdict(set)}
    
    for gene_id, go_terms in go_data.items():
        for go_id, go_info in go_terms.items():
            go_type = go_info['type']
            if go_type in go_by_type:
                go_by_type[go_type][go_id].add(gene_id)
    
    return go_by_type

def perform_enrichment_analysis(study_genes, population_genes, go_by_type, go_type):
    """Perform Fisher's exact test for GO enrichment"""
    results = []
    
    study_set = set(study_genes)
    population_set = set(population_genes)
    
    # Only consider genes that are in both datasets
    common_genes = study_set.intersection(population_set)
    study_common = len(common_genes)
    population_common = len(population_set)
    
    print(f"GO Type {go_type}: {study_common} study genes, {population_common} population genes")
    
    for go_id, annotated_genes in go_by_type[go_type].items():
        # Genes annotated with this GO term in the population
        annotated_pop = annotated_genes.intersection(population_set)
        # Genes annotated with this GO term in the study set
        annotated_study = annotated_genes.intersection(study_set)
        
        if len(annotated_pop) == 0:
            continue
        
        # 2x2 contingency table
        # |  GO+  |  GO-  |
        # | a  b  | Study |
        # | c  d  | Pop   |
        
        a = len(annotated_study)  # Study genes with GO term
        b = study_common - a      # Study genes without GO term
        c = len(annotated_pop) - a  # Population genes with GO term (excluding study)
        d = population_common - len(annotated_pop) - b  # Population genes without GO term (excluding study)
        
        if a == 0:  # No enrichment if no study genes have the term
            continue
            
        # Fisher's exact test
        odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
        
        # Calculate fold enrichment
        expected = (len(annotated_pop) / population_common) * study_common
        fold_enrichment = a / expected if expected > 0 else float('inf')
        
        results.append({
            'GO_ID': go_id,
            'GO_Type': go_type,
            'Study_Count': a,
            'Study_Total': study_common,
            'Population_Count': len(annotated_pop),
            'Population_Total': population_common,
            'Fold_Enrichment': fold_enrichment,
            'P_Value': p_value,
            'Odds_Ratio': odds_ratio
        })
    
    return results

def main():
    parser = argparse.ArgumentParser(description='GO Enrichment Analysis')
    parser.add_argument('--study', required=True, help='Study genes file (gene_mapping.tsv)')
    parser.add_argument('--population', required=True, help='Population genes file (full_genome_go_map.tsv)')
    parser.add_argument('--output', default='go_enrichment_results', help='Output prefix')
    parser.add_argument('--alpha', type=float, default=0.05, help='FDR significance threshold')
    
    args = parser.parse_args()
    
    print("Loading study genes...")
    study_go_data = load_go_data(args.study)
    study_genes = list(study_go_data.keys())
    
    print("Loading population genes...")
    population_go_data = load_go_data(args.population)
    population_genes = list(population_go_data.keys())
    
    print(f"Study genes: {len(study_genes)}")
    print(f"Population genes: {len(population_genes)}")
    
    # Separate GO terms by type for both datasets
    print("Separating GO terms by type...")
    study_go_by_type = separate_go_by_type(study_go_data)
    population_go_by_type = separate_go_by_type(population_go_data)
    
    # Combine GO terms from both datasets for analysis
    combined_go_by_type = {'C': defaultdict(set), 'F': defaultdict(set), 'P': defaultdict(set)}
    for go_type in ['C', 'F', 'P']:
        for go_id in set(list(study_go_by_type[go_type].keys()) + list(population_go_by_type[go_type].keys())):
            combined_go_by_type[go_type][go_id] = (
                study_go_by_type[go_type][go_id].union(population_go_by_type[go_type][go_id])
            )
    
    # Perform enrichment analysis for each GO type
    all_results = []
    go_type_names = {'C': 'Cellular_Component', 'F': 'Molecular_Function', 'P': 'Biological_Process'}
    
    for go_type in ['C', 'F', 'P']:
        print(f"\nPerforming enrichment analysis for {go_type} ({go_type_names[go_type]})...")
        results = perform_enrichment_analysis(study_genes, population_genes, combined_go_by_type, go_type)
        
        if results:
            # Add GO term names
            for result in results:
                go_id = result['GO_ID']
                # Try to get term name from study data first, then population
                term_name = 'Unknown'
                for gene_id, go_terms in study_go_data.items():
                    if go_id in go_terms:
                        term_name = go_terms[go_id]['name']
                        break
                if term_name == 'Unknown':
                    for gene_id, go_terms in population_go_data.items():
                        if go_id in go_terms:
                            term_name = go_terms[go_id]['name']
                            break
                result['GO_Term'] = term_name
            
            all_results.extend(results)
            
            # Save type-specific results
            df_type = pd.DataFrame(results)
            if not df_type.empty:
                # Apply FDR correction
                df_type['FDR_P_Value'] = false_discovery_control(df_type['P_Value'].values)
                df_type = df_type.sort_values('P_Value')
                
                output_file = f"{args.output}_{go_type_names[go_type].lower()}.tsv"
                df_type.to_csv(output_file, sep='\t', index=False)
                print(f"Saved {len(df_type)} results to {output_file}")
                
                # Print significant results
                significant = df_type[df_type['FDR_P_Value'] < args.alpha]
                if not significant.empty:
                    print(f"Significant terms (FDR < {args.alpha}): {len(significant)}")
                    for _, row in significant.head(10).iterrows():
                        print(f"  {row['GO_ID']}: {row['GO_Term']} (p={row['P_Value']:.2e}, FDR={row['FDR_P_Value']:.2e})")
    
    # Save combined results
    if all_results:
        df_all = pd.DataFrame(all_results)
        df_all['FDR_P_Value'] = false_discovery_control(df_all['P_Value'].values)
        df_all = df_all.sort_values('P_Value')
        
        output_file = f"{args.output}_combined.tsv"
        df_all.to_csv(output_file, sep='\t', index=False)
        print(f"\nSaved {len(df_all)} total results to {output_file}")
        
        # Summary statistics
        for go_type in ['C', 'F', 'P']:
            type_results = df_all[df_all['GO_Type'] == go_type]
            significant = type_results[type_results['FDR_P_Value'] < args.alpha]
            print(f"{go_type_names[go_type]}: {len(type_results)} terms tested, {len(significant)} significant")

if __name__ == "__main__":
    main()