from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def lr_distribution_workflow(config_file: str, gwf):
    """
    Create workflow for likelihood ratio distribution analysis.
    
    Args:
        config_file: Path to configuration YAML file
        gwf: GWF workflow object
    """
    
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
    SWEEPFINDER_RESULT: str = CONFIG['sweepfinder2_result_file']
    
    # --------------------------------------------------
    #          Likelihood Ratio Distribution Analysis
    # --------------------------------------------------
    
    # Create likelihood ratio distribution analysis target
    gwf.target_from_template(
        name = f"{SPECIES_ID}_lr_distribution_analysis",
        template = lr_distribution_analysis_template(
            work_path = WORK_DIR,
            script_path = SCRIPTS_PATH,
            log_path = LOG_DIR,
            spid = SPECIES_ID,
            input_file = SWEEPFINDER_RESULT,
            output_dir = f"{OUTPUT_DIR}/{SPECIES_ID}_lr_distribution"
        )
    )
    
    return gwf
