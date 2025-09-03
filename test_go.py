import requests
import pandas as pd

def get_go_terms(uniprot_id):
      """Get GO terms for a UniProt accession"""
      # Clean the ID (remove tr|, sp| prefixes)
      clean_id = uniprot_id.split('|')[1] if '|' in uniprot_id else uniprot_id

      url = f"https://rest.uniprot.org/uniprotkb/{clean_id}.json"

      try:
          response = requests.get(url)
          if response.status_code == 200:
              data = response.json()

              # Extract GO terms
              go_terms = []
              if 'dbReferences' in data:
                  for ref in data['dbReferences']:
                      if ref['type'] == 'GO':
                          go_id = ref['id']
                          go_term = ref['properties']['term']
                          go_evidence = ref['properties']['evidence']
                          go_terms.append({
                              'GO_ID': go_id,
                              'GO_Term': go_term,
                              'Evidence': go_evidence
                          })

              return go_terms
      except:
          return []

# Example usage
uniprot_id = "A0A123456"  # From your Diamond results
go_terms = get_go_terms(uniprot_id)
print(go_terms)
