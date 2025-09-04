from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def region_gene_workflow(config_file: str, gwf):

    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    CONFIG = yaml.safe_load(open(config_file))
    ACCOUNT: str = CONFIG['account']
    SPECIES_ID: str = CONFIG['species_id']
    WORK_DIR: str = CONFIG['working_directory_path']
    OUTPUT_DIR: str = CONFIG['output_directory_path']
    LOG_DIR: str = CONFIG['log_directory_path']
    SCRIPTS_PATH: str = CONFIG['scripts_path']
    
    # Region gene analysis parameters
    REGIONS_FILE: str = CONFIG['regions_file']
    GENOME_FILE: str = CONFIG['genome_file']
    GFF_FILE: str = CONFIG['gff_file']
    MIN_LENGTH: int = CONFIG.get('min_region_length', 5000)
    FLANK_SIZE: int = CONFIG.get('flank_size', 10000)
    
    # BLAST parameters
    BLAST_DB: str = CONFIG['blast_db']
    EVALUE: float = CONFIG.get('blast_evalue', 1e-5)
    
    # Whole genome processing option
    PROCESS_WHOLE_GENOME: bool = CONFIG.get('process_whole_genome', False)
    GENE_CHUNK_SIZE: int = CONFIG.get('gene_chunk_size', 1000)
    
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    # --------------------------------------------------
    #          Extract genes from all regions
    # --------------------------------------------------
    gwf.target_from_template(
        name = f"{SPECIES_ID}_extract_region_genes",
        template = extract_region_genes(
            work_path = WORK_DIR, 
            script_path = SCRIPTS_PATH, 
            log_path = LOG_DIR, 
            spid = SPECIES_ID,
            regions_file = REGIONS_FILE,
            genome_file = GENOME_FILE,
            gff_file = GFF_FILE,
            min_length = MIN_LENGTH,
            flank_size = FLANK_SIZE) 
    )

    # --------------------------------------------------
    #          Diamond gene annotation (parallel by region)
    # --------------------------------------------------
    
    # Get unique regions to create parallel Diamond jobs
    def get_unique_regions(regions_file, min_length):
        """Parse regions file and return unique regions meeting length threshold"""
        unique_regions = {}
        
        try:
            with open(regions_file, 'r') as f:
                header = f.readline().strip().split('\t')
                region_id_col = header.index('region_id')
                start_col = header.index('start') 
                end_col = header.index('end')
                
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) > max(region_id_col, start_col, end_col):
                        region_id = parts[region_id_col]
                        start = int(parts[start_col])
                        end = int(parts[end_col])
                        length = end - start + 1
                        
                        # Only include regions meeting length threshold
                        if length >= min_length and region_id not in unique_regions:
                            unique_regions[region_id] = {'start': start, 'end': end, 'length': length}
        
        except Exception as e:
            print(f"Warning: Could not parse regions file: {e}")
            return {}
        
        return unique_regions
    
    # Create parallel Diamond jobs for each region
    unique_regions = get_unique_regions(REGIONS_FILE, MIN_LENGTH)
    
    if unique_regions:
        for region_id in unique_regions:
            gwf.target_from_template(
                name = f"{SPECIES_ID}_diamond_{region_id}",
                template = blast_gene_annotation(
                    work_path = WORK_DIR,
                    script_path = SCRIPTS_PATH,
                    log_path = LOG_DIR,
                    spid = SPECIES_ID,
                    blast_db = BLAST_DB,
                    region_id = region_id,
                    evalue = EVALUE)
            )
    else:
        # Fallback to single job for all regions
        gwf.target_from_template(
            name = f"{SPECIES_ID}_diamond_gene_annotation",
            template = blast_gene_annotation(
                work_path = WORK_DIR,
                script_path = SCRIPTS_PATH,
                log_path = LOG_DIR,
                spid = SPECIES_ID,
                blast_db = BLAST_DB,
                evalue = EVALUE)
        )
    
    # --------------------------------------------------
    #          Whole genome gene processing (optional)
    # --------------------------------------------------
    if PROCESS_WHOLE_GENOME:
        
        # Extract all genes from genome into single FASTA file
        gwf.target_from_template(
            name = f"{SPECIES_ID}_extract_all_genes",
            template = extract_all_genes(
                work_path = WORK_DIR, 
                script_path = SCRIPTS_PATH, 
                log_path = LOG_DIR, 
                spid = SPECIES_ID,
                genome_file = GENOME_FILE,
                gff_file = GFF_FILE) 
        )

        # Run Diamond annotation on all genes
        gwf.target_from_template(
            name = f"{SPECIES_ID}_diamond_all_genes",
            template = blast_all_genes_annotation(
                work_path = WORK_DIR,
                script_path = SCRIPTS_PATH,
                log_path = LOG_DIR,
                spid = SPECIES_ID,
                blast_db = BLAST_DB,
                evalue = EVALUE)
        )

    return gwf