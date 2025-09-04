import requests
import pandas as pd
import argparse
import sys
import time
from typing import List, Dict

def read_accession_file(file_path: str) -> List[str]:
    """Read UniProt accession IDs from a text file"""
    accessions = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    accessions.append(line)
    except FileNotFoundError:
        print(f"Error: File {file_path} not found", file=sys.stderr)
        sys.exit(1)
    
    return accessions

def submit_mapping_job(accessions: List[str], from_db: str = "UniProtKB_AC-ID", to_db: str = "UniProtKB") -> str:
    """Submit ID mapping job to UniProt"""
    url = "https://rest.uniprot.org/idmapping/run"
    
    # Join accessions with commas (like the curl command)
    ids_string = ','.join(accessions)
    
    data = {
        'from': from_db,
        'to': to_db,
        'ids': ids_string
    }
    
    print(f"Submitting mapping job for {len(accessions)} IDs...", file=sys.stderr)
    
    try:
        response = requests.post(url, data=data)
        response.raise_for_status()
        
        result = response.json()
        job_id = result['jobId']
        print(f"Job submitted successfully. Job ID: {job_id}", file=sys.stderr)
        return job_id
        
    except requests.exceptions.RequestException as e:
        print(f"Error submitting job: {e}", file=sys.stderr)
        sys.exit(1)

def check_job_status(job_id: str) -> str:
    """Check the status of a mapping job"""
    url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    max_retries = 60  # Maximum number of status checks (5 minutes)
    retry_count = 0
    
    while retry_count < max_retries:
        try:
            response = requests.get(url)
            response.raise_for_status()
            
            status_data = response.json()
            print(f"Status check {retry_count + 1}: {status_data}", file=sys.stderr)
            
            if 'jobStatus' in status_data:
                status = status_data['jobStatus']
                if status in ['NEW', 'RUNNING']:
                    print(f"Job status: {status}, waiting...", file=sys.stderr)
                    time.sleep(5)
                    retry_count += 1
                    continue
                elif status == 'FINISHED':
                    print("Job completed successfully!", file=sys.stderr)
                    return 'FINISHED'
                else:
                    print(f"Job failed with status: {status}", file=sys.stderr)
                    return status
            else:
                # No jobStatus field means job is completed
                if 'results' in status_data or 'failedIds' in status_data:
                    print("Job completed successfully!", file=sys.stderr)
                    return 'FINISHED'
                else:
                    print(f"Unknown status response: {status_data}", file=sys.stderr)
                    retry_count += 1
                    time.sleep(5)
                    continue
                
        except requests.exceptions.RequestException as e:
            print(f"Error checking job status: {e}", file=sys.stderr)
            retry_count += 1
            time.sleep(5)
            continue
    
    print("Timeout waiting for job completion", file=sys.stderr)
    return 'TIMEOUT'

def get_mapping_results(job_id: str) -> Dict:
    """Retrieve mapping results"""
    url = f"https://rest.uniprot.org/idmapping/stream/{job_id}"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        
        results = response.json()
        print(f"Retrieved mapping results", file=sys.stderr)
        return results
        
    except requests.exceptions.RequestException as e:
        print(f"Error retrieving results: {e}", file=sys.stderr)
        return {}

def fetch_protein_details(accession: str, debug_log_path: str = None) -> tuple:
    """Fetch detailed protein information from UniProt"""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            
            if debug_log_path:
                with open(debug_log_path, 'a') as f:
                    f.write(f"\n--- Fetching details for {accession} ---\n")
                    f.write(f"Keys available: {list(data.keys())}\n")
                    
                    # Explore uniProtKBCrossReferences structure
                    if 'uniProtKBCrossReferences' in data:
                        f.write(f"Found {len(data['uniProtKBCrossReferences'])} cross references\n")
                        
                        # Show first few cross references
                        for i, ref in enumerate(data['uniProtKBCrossReferences'][:5]):
                            f.write(f"Cross-ref {i+1}: {ref}\n")
                        
                        # Show GO-specific references
                        go_refs = [ref for ref in data['uniProtKBCrossReferences'] if ref.get('database') == 'GO']
                        f.write(f"GO references found: {len(go_refs)}\n")
                        for go_ref in go_refs[:3]:
                            f.write(f"GO example: {go_ref}\n")
                    else:
                        f.write("No uniProtKBCrossReferences found\n")
            
            # Extract information
            entry_name = data.get('uniProtkbId', 'NA')
            protein_name = 'NA'
            gene_names = 'NA'
            organism = 'NA'
            
            # Extract protein name
            if 'proteinDescription' in data:
                desc = data['proteinDescription']
                if 'recommendedName' in desc and 'fullName' in desc['recommendedName']:
                    protein_name = desc['recommendedName']['fullName'].get('value', 'NA')
                elif 'submissionNames' in desc and desc['submissionNames']:
                    if 'fullName' in desc['submissionNames'][0]:
                        protein_name = desc['submissionNames'][0]['fullName'].get('value', 'NA')
            
            # Extract gene names
            if 'genes' in data and data['genes']:
                gene_list = []
                for gene in data['genes']:
                    if 'geneName' in gene:
                        gene_list.append(gene['geneName']['value'])
                    if 'synonyms' in gene:
                        gene_list.extend([syn['value'] for syn in gene['synonyms']])
                if gene_list:
                    gene_names = ';'.join(gene_list)
            
            # Extract organism
            if 'organism' in data:
                organism = data['organism'].get('scientificName', 'NA')
            
            # Extract GO terms from uniProtKBCrossReferences
            go_ids = []
            go_terms = []
            go_evidence = []
            
            if 'uniProtKBCrossReferences' in data:
                for ref in data['uniProtKBCrossReferences']:
                    if ref.get('database') == 'GO':
                        go_id = ref.get('id', 'NA')
                        go_term = 'NA'
                        evidence = 'NA'
                        
                        # Extract GO term description and evidence
                        if 'properties' in ref:
                            for prop in ref['properties']:
                                if prop.get('key') == 'GoTerm':
                                    go_term = prop.get('value', 'NA')
                                elif prop.get('key') == 'GoEvidenceType':
                                    evidence = prop.get('value', 'NA')
                        
                        go_ids.append(go_id)
                        go_terms.append(go_term)
                        go_evidence.append(evidence)
            
            # Format GO information (semicolon-separated)
            go_ids_str = ';'.join(go_ids) if go_ids else 'NA'
            go_terms_str = ';'.join(go_terms) if go_terms else 'NA'
            go_evidence_str = ';'.join(go_evidence) if go_evidence else 'NA'
            
            if debug_log_path:
                with open(debug_log_path, 'a') as f:
                    f.write(f"Details extracted: Name={entry_name}, Protein={protein_name}, Genes={gene_names}, Organism={organism}\n")
                    f.write(f"GO terms: {len(go_ids)} found - IDs={go_ids_str[:100]}...\n")
            
            return entry_name, protein_name, gene_names, organism, go_ids_str, go_terms_str, go_evidence_str
        
        else:
            if debug_log_path:
                with open(debug_log_path, 'a') as f:
                    f.write(f"Failed to fetch {accession}: HTTP {response.status_code}\n")
            return 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'
            
    except Exception as e:
        if debug_log_path:
            with open(debug_log_path, 'a') as f:
                f.write(f"Error fetching {accession}: {e}\n")
        return 'NA', 'NA', 'NA', 'NA'

def process_results_to_dataframe(results: Dict, original_ids: List[str], debug_log_path: str = None) -> pd.DataFrame:
    """Convert mapping results to DataFrame preserving all original IDs including duplicates"""
    
    # Create mapping dictionary from API results
    mapping_dict = {}
    
    if 'results' in results:
        # Debug: save sample of results to file (if debug log specified)
        if debug_log_path:
            with open(debug_log_path, 'a') as f:
                f.write(f"Total API results: {len(results['results'])}\n")
                f.write(f"Sample result structure:\n")
                if results['results']:
                    f.write(f"First result: {results['results'][0]}\n")
                    f.write(f"Type of 'to' field: {type(results['results'][0]['to'])}\n")
        
        # Build mapping dictionary from API results
        for result in results['results']:
            from_id = result['from']
            to_entry = result['to']
            
            # Debug: log what type we got (if debug log specified)
            if debug_log_path:
                with open(debug_log_path, 'a') as f:
                    f.write(f"\nProcessing API result {from_id}: type={type(to_entry)}, value={str(to_entry)[:200]}...\n")
            
            # Handle different types of 'to' entries
            if isinstance(to_entry, str):
                # Simple string result (e.g., just accession ID)
                # Fetch detailed information for this accession
                entry_name, protein_name, gene_names, organism, go_ids, go_terms, go_evidence = fetch_protein_details(to_entry, debug_log_path)
                
                mapping_dict[from_id] = {
                    'Entry': to_entry,
                    'Entry_Name': entry_name,
                    'Protein_Name': protein_name,
                    'Gene_Names': gene_names,
                    'Organism': organism,
                    'GO_IDs': go_ids,
                    'GO_Terms': go_terms,
                    'GO_Evidence': go_evidence,
                    'Status': 'Success'
                }
            elif isinstance(to_entry, dict):
                # Full entry object - save structure for debugging (if debug log specified)
                if debug_log_path:
                    debug_info = f"Processing entry keys: {list(to_entry.keys())}\n"
                    debug_info += f"Full entry sample: {str(to_entry)[:500]}...\n"
                    
                    with open(debug_log_path, 'a') as debug_file:
                        debug_file.write(f"\n=== Processing {from_id} ===\n")
                        debug_file.write(debug_info)
                
                uniprot_id = to_entry.get('primaryAccession', 'NA')
                entry_name = to_entry.get('uniProtkbId', 'NA')
                protein_name = 'NA'
                gene_names = 'NA'
                organism = 'NA'
                
                # Extract protein name - try multiple possible structures
                if 'proteinDescription' in to_entry:
                    desc = to_entry['proteinDescription']
                    if isinstance(desc, dict):
                        if 'recommendedName' in desc:
                            rec_name = desc['recommendedName']
                            if isinstance(rec_name, dict) and 'fullName' in rec_name:
                                if isinstance(rec_name['fullName'], dict):
                                    protein_name = rec_name['fullName'].get('value', 'NA')
                                else:
                                    protein_name = str(rec_name['fullName'])
                        elif 'submissionNames' in desc and desc['submissionNames']:
                            sub_name = desc['submissionNames'][0]
                            if isinstance(sub_name, dict) and 'fullName' in sub_name:
                                if isinstance(sub_name['fullName'], dict):
                                    protein_name = sub_name['fullName'].get('value', 'NA')
                                else:
                                    protein_name = str(sub_name['fullName'])
                
                # Extract gene names - handle different structures
                if 'genes' in to_entry and to_entry['genes']:
                    gene_list = []
                    for gene in to_entry['genes']:
                        if isinstance(gene, dict):
                            if 'geneName' in gene:
                                if isinstance(gene['geneName'], dict):
                                    gene_list.append(gene['geneName'].get('value', ''))
                                else:
                                    gene_list.append(str(gene['geneName']))
                            if 'synonyms' in gene and gene['synonyms']:
                                for syn in gene['synonyms']:
                                    if isinstance(syn, dict):
                                        gene_list.append(syn.get('value', ''))
                                    else:
                                        gene_list.append(str(syn))
                    # Filter out empty strings and join
                    gene_list = [g for g in gene_list if g]
                    if gene_list:
                        gene_names = ';'.join(gene_list)
                
                # Extract organism - handle different structures
                if 'organism' in to_entry:
                    org = to_entry['organism']
                    if isinstance(org, dict):
                        organism = org.get('scientificName', 'NA')
                    else:
                        organism = str(org)
                
                # Save extraction results to debug file (if debug log specified)
                if debug_log_path:
                    with open(debug_log_path, 'a') as debug_file:
                        debug_file.write(f"Extracted: ID={uniprot_id}, Name={entry_name}, Protein={protein_name}, Genes={gene_names}, Organism={organism}\n")
                
                mapping_dict[from_id] = {
                    'Entry': uniprot_id,
                    'Entry_Name': entry_name,
                    'Protein_Name': protein_name,
                    'Gene_Names': gene_names,
                    'Organism': organism,
                    'GO_IDs': 'NA',
                    'GO_Terms': 'NA',
                    'GO_Evidence': 'NA',
                    'Status': 'Success'
                }
            else:
                # Unknown format
                mapping_dict[from_id] = {
                    'Entry': str(to_entry),
                    'Entry_Name': 'NA',
                    'Protein_Name': 'NA',
                    'Gene_Names': 'NA',
                    'Organism': 'NA',
                    'GO_IDs': 'NA',
                    'GO_Terms': 'NA',
                    'GO_Evidence': 'NA',
                    'Status': 'Success'
                }
    
    # Now process ALL original IDs (preserving duplicates)
    mapping_data = []
    mapped_count = 0
    not_found_count = 0
    
    for original_id in original_ids:
        clean_id = original_id.split('|')[1] if '|' in original_id else original_id
        
        if clean_id in mapping_dict:
            # Found mapping - use the data
            mapping_info = mapping_dict[clean_id]
            mapping_data.append({
                'Query': clean_id,
                **mapping_info
            })
            mapped_count += 1
        else:
            # Not found - add NA entry
            mapping_data.append({
                'Query': clean_id,
                'Entry': 'NA',
                'Entry_Name': 'NA',
                'Protein_Name': 'NA',
                'Gene_Names': 'NA',
                'Organism': 'NA',
                'GO_IDs': 'NA',
                'GO_Terms': 'NA',
                'GO_Evidence': 'NA',
                'Status': 'Not_Found'
            })
            not_found_count += 1
    
    # Log mapping summary
    if debug_log_path:
        with open(debug_log_path, 'a') as f:
            f.write(f"\n=== MAPPING SUMMARY ===\n")
            f.write(f"Input IDs (total): {len(original_ids)}\n")
            f.write(f"Input IDs (unique): {len(set(original_ids))}\n")
            f.write(f"Successfully mapped: {mapped_count}\n")
            f.write(f"Not found: {not_found_count}\n")
            f.write(f"Output rows: {len(mapping_data)}\n")
            f.write(f"Input == Output: {len(original_ids) == len(mapping_data)}\n")
    
    return pd.DataFrame(mapping_data)

def main():
    parser = argparse.ArgumentParser(
        description='UniProt ID batch mapping using REST API',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python uniprot_batch_mapper.py -i accessions.txt -o results.tsv
  python uniprot_batch_mapper.py -i ids.txt -f UniProtKB_AC-ID -t Gene_Name -o gene_mapping.tsv
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input text file with UniProt accession IDs (one per line)')
    parser.add_argument('-o', '--output', default='uniprot_mapping.tsv',
                       help='Output TSV file (default: uniprot_mapping.tsv)')
    parser.add_argument('-f', '--from_db', default='UniProtKB_AC-ID',
                       help='Source database/format (default: UniProtKB_AC-ID)')
    parser.add_argument('-t', '--to_db', default='UniProtKB',
                       help='Target database/format (default: UniProtKB)')
    parser.add_argument('-d', '--debug_log',
                       help='Path to save debug log file (optional)')
    
    args = parser.parse_args()
    
    # Read accession IDs from file
    accessions = read_accession_file(args.input)
    
    if not accessions:
        print("No accession IDs found in input file", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(accessions)} accession IDs", file=sys.stderr)
    
    # Optional debug analysis
    if args.debug_log:
        unique_accessions = list(set(accessions))
        duplicates = len(accessions) - len(unique_accessions)
        
        with open(args.debug_log, 'w') as f:
            f.write(f"=== INPUT ANALYSIS ===\n")
            f.write(f"Total input IDs: {len(accessions)}\n")
            f.write(f"Unique input IDs: {len(unique_accessions)}\n")
            f.write(f"Duplicate IDs: {duplicates}\n")
            f.write(f"First 10 input IDs: {accessions[:10]}\n")
            f.write(f"Last 10 input IDs: {accessions[-10:]}\n\n")
    
    # Submit mapping job
    job_id = submit_mapping_job(accessions, args.from_db, args.to_db)
    
    # Wait for job completion
    status = check_job_status(job_id)
    
    if status != 'FINISHED':
        print(f"Job failed with status: {status}", file=sys.stderr)
        sys.exit(1)
    
    # Get results
    results = get_mapping_results(job_id)
    
    if not results:
        print("No results retrieved", file=sys.stderr)
        sys.exit(1)
    
    # Convert to DataFrame
    df = process_results_to_dataframe(results, accessions, args.debug_log)
    
    # Save to file
    df.to_csv(args.output, sep='\t', index=False)
    
    # Print summary
    success_count = len(df[df['Status'] == 'Success'])
    not_found_count = len(df[df['Status'] == 'Not_Found'])
    
    print(f"Results saved to: {args.output}", file=sys.stderr)
    print(f"Successfully mapped: {success_count}", file=sys.stderr)
    print(f"Not found: {not_found_count}", file=sys.stderr)
    print(f"Total processed: {len(df)}", file=sys.stderr)

if __name__ == "__main__":
    main()
