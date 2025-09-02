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

    return gwf