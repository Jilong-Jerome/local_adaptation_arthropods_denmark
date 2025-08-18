from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def sweepfinder2_workflow(config_file: str, gwf):

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
    POPS: str = CONFIG['pop_order_in_vcf']
    CHRS: str = CONFIG['reference_fai_index']
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    # --------------------------------------------------
    #          Parse input for sweepfinder2
    # --------------------------------------------------
    for maf in [0.05]:
        gwf.target_from_template(
            name = f"{SPECIES_ID}_sweepfinder2_parse_input_maf{maf}_jilong",
            template = parse_input_sweepfinder2_jilong(
                work_path = WORK_DIR, 
                script_path =  SCRIPTS_PATH, 
                log_path =  LOG_DIR, 
                spid =  SPECIES_ID, 
                vcf = INPUT_VCF, 
                maf = maf) 
        )
        pops = open(POPS,"r")
        pop_col = 3
        pop_size = 100
        pop_logs = {}
        for pop in pops:
            pop = pop.strip("\n")
            pop_id = pop.replace("-","_")
            chroms = open(CHRS,"r")
            for line in chroms:
                chrom = line.strip("\n").split("\t")[0]
                chrom_len = int(line.strip("\n").split("\t")[1])
                chrom_id = chrom.replace("|","_")
                gwf.target_from_template(
                    name = f"{SPECIES_ID}_split_{pop_id}_{chrom_id}_maf{maf}",
                    template = split_sweefinder_input(
                        work_path = WORK_DIR,
                        script_path = SCRIPTS_PATH,
                        log_path = LOG_DIR,
                        spid = SPECIES_ID,
                        maf = maf,
                        chrom = chrom,
                        chrom_id = chrom_id,
                        pop = pop,
                        pop_col = pop_col,
                        pop_size = pop_size) 
                )
                gwf.target_from_template(
                    name = f"{SPECIES_ID}_specturm_{pop_id}_{chrom_id}_maf{maf}",
                    template = sweepfinder_chr_pop_spectrum(
                            work_path = WORK_DIR,
                            script_path = SCRIPTS_PATH,
                            log_path = LOG_DIR,
                            spid = SPECIES_ID,
                            maf = maf,
                            chrom = chrom,
                            chrom_id = chrom_id,
                            pop = pop)
                    )
                for grid in [1000]:
                    for option in [1,2]:
                        if (grid,option) not in pop_logs:
                            pop_logs[(grid,option)] = []
                        grid_size = grid
                        gwf.target_from_template(
                            name = f"{SPECIES_ID}_sweepfinder2_{pop_id}_{chrom_id}_maf{maf}_w{grid_size}_o{option}",
                            template = sweepfinder_chr_pop_scan(
                                work_path = WORK_DIR,
                                script_path = SCRIPTS_PATH,
                                log_path = LOG_DIR,
                                spid = SPECIES_ID,
                                maf = maf,
                                chrom = chrom,
                                chrom_id = chrom_id,
                                pop = pop,
                                grid_size = grid_size,
                                option = option) 
                        )
                        pop_logs[(grid,option)].append(f"{LOG_DIR}/{SPECIES_ID}/{SPECIES_ID}_{chrom_id}_{pop}_{maf}_w{grid_size}_o{option}_sweepfinder2.DONE")
            pop_col = pop_col + 1
        # Concatenate sweepfinder results
        for grid in [1000]:
            for option in [1,2]:
                gwf.target_from_template(
                    name = f"{SPECIES_ID}_concat_sweepfinder2_maf{maf}_w{grid}_o{option}",
                    template = sweepfinder_concat_res(
                        work_path = WORK_DIR,
                        script_path = SCRIPTS_PATH,
                        log_path = LOG_DIR,
                        spid = SPECIES_ID,
                        maf = maf,
                        grid_size = grid,
                        option = option,
                        chr_pos_logs = pop_logs[(grid,option)])
                )
                gwf.target_from_template(
                    name = f"{SPECIES_ID}_manhattan_maf{maf}_w{grid}_o{option}",
                    template = sweepfinder_viz_res(
                        work_path = WORK_DIR,
                        script_path = SCRIPTS_PATH,
                        log_path = LOG_DIR,
                        spid = SPECIES_ID,
                        maf = maf,
                        grid_size = grid,
                        option = option)
                )
                    

    return gwf
