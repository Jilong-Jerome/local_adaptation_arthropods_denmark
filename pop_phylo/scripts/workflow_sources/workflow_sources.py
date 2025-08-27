from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def pop_phylo_workflow(config_file: str, gwf):

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
    INPUT_VCF: str = CONFIG['input_vcf']
    
    # Phylogeny-specific parameters
    BOOTSTRAP_REPS: int = CONFIG.get('bootstrap_replicates', 100)
    SEED: int = CONFIG.get('random_seed', 42)
    OUTPUT_NAME: str = CONFIG.get('output_name', 'phylo_sequences')
    SCAFFOLDS: str = CONFIG.get('scaffolds', '')
    REGIONS: str = CONFIG.get('regions', '')
    
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    # --------------------------------------------------
    #          Bootstrap sequence sampling
    # --------------------------------------------------
    gwf.target_from_template(
        name = f"{SPECIES_ID}_bootstrap_sampling",
        template = bootstrap_sampling(
            work_path = WORK_DIR, 
            script_path = SCRIPTS_PATH, 
            log_path = LOG_DIR, 
            spid = SPECIES_ID, 
            vcf = INPUT_VCF, 
            bootstrap_reps = BOOTSTRAP_REPS,
            seed = SEED,
            output_name = OUTPUT_NAME,
            scaffolds = SCAFFOLDS,
            regions = REGIONS) 
    )

    return gwf