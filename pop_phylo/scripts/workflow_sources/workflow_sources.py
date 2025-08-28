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
    FST_FILE: str = CONFIG.get('fst_file', '')
    FST_METHOD: str = CONFIG.get('fst_method', 'average')
    KEEP_ORIGINAL_NAMES: bool = CONFIG.get('keep_original_names', True)
    IQTREE_MODEL: str = CONFIG.get('iqtree_model', 'MFP+G')
    IQTREE_THREADS: int = CONFIG.get('iqtree_threads', 4)
    
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

    # --------------------------------------------------
    #          FST-based phylogeny generation
    # --------------------------------------------------
    if FST_FILE:  # Only add this step if FST file is provided
        gwf.target_from_template(
            name = f"{SPECIES_ID}_fst_phylogeny",
            template = fst_phylo_generation(
                work_path = WORK_DIR, 
                script_path = SCRIPTS_PATH, 
                log_path = LOG_DIR, 
                spid = SPECIES_ID, 
                fst_file = FST_FILE,
                output_name = OUTPUT_NAME,
                method = FST_METHOD) 
        )

    # --------------------------------------------------
    #          IQ-TREE input preparation (parallel)
    # --------------------------------------------------
    if FST_FILE:  # Only add these steps if FST file is provided
        for rep in range(1, BOOTSTRAP_REPS + 1):
            gwf.target_from_template(
                name = f"{SPECIES_ID}_iqtree_prep_rep{rep:03d}",
                template = prepare_iqtree_single_rep(
                    work_path = WORK_DIR, 
                    script_path = SCRIPTS_PATH, 
                    log_path = LOG_DIR, 
                    spid = SPECIES_ID, 
                    replicate = rep,
                    output_name = OUTPUT_NAME,
                    keep_original_names = KEEP_ORIGINAL_NAMES) 
            )

    # --------------------------------------------------
    #          IQ-TREE3 phylogenetic inference (parallel)
    # --------------------------------------------------
    if FST_FILE:  # Only add these steps if FST file is provided
        for rep in range(1, BOOTSTRAP_REPS + 1):
            gwf.target_from_template(
                name = f"{SPECIES_ID}_iqtree_rep{rep:03d}",
                template = run_iqtree3_single_rep(
                    work_path = WORK_DIR, 
                    script_path = SCRIPTS_PATH, 
                    log_path = LOG_DIR, 
                    spid = SPECIES_ID, 
                    replicate = rep,
                    model = IQTREE_MODEL,
                    threads = IQTREE_THREADS) 
            )

    return gwf