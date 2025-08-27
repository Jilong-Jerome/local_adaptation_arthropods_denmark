#!/bin/env python3
from gwf import AnonymousTarget
import os, glob


def bootstrap_sampling(work_path: str, script_path: str, log_path: str, spid: str, vcf: str, 
                      bootstrap_reps: int, seed: int, output_name: str, scaffolds: str = "", 
                      regions: str = ""):
    """Generate bootstrap sequence samples from VCF for phylogenetic analysis."""
    inputs = {"vcf": vcf}
    
    scaffold_args = f"-s {scaffolds}" if scaffolds else ""
    region_args = f"--regions {regions}" if regions else ""
    seed_arg = f"--seed {seed}" if seed else ""
    
    outputs = {
        "log": f"{log_path}/{spid}/{spid}_bootstrap_sampling.DONE"
    }
    
    options = {
        'cores': 1,
        'memory': '8g',
        'walltime': '12:00:00',
        'account': 'EcoGenetics'
    }
    
    spec = f"""
    # Setting conda environments
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate python_phylo
    
    # Writting job information to standard output
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    # Setting working directory
    mkdir -p {work_path}/{spid}/bootstrap_sequences
    mkdir -p {log_path}/{spid}
    
    cd {work_path}/{spid}/bootstrap_sequences
    
    # Run bootstrap sampling using our existing script
    python {script_path}/bootstrap_phylo_sampler.py \\
        {inputs["vcf"]} \\
        --outdir {work_path}/{spid}/bootstrap_sequences \\
        -o {output_name} \\
        -r {bootstrap_reps} \\
        {scaffold_args} \\
        {region_args} \\
        {seed_arg}
    
    echo "FINISH: $(date)"
    jobinfo $SLURM_JOBID
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
