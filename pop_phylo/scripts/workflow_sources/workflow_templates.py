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
        'walltime': '48:00:00',
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


def fst_phylo_generation(work_path: str, script_path: str, log_path: str, spid: str, 
                        fst_file: str, output_name: str, method: str = "average"):
    """Generate FST-based phylogeny from FST distance matrix."""
    inputs = {"fst_file": fst_file}
    
    outputs = {
        "newick": f"{work_path}/{spid}/{spid}_fst_phylo.newick",
        "log": f"{log_path}/{spid}/{spid}_fst_phylo.DONE"
    }
    
    options = {
        'cores': 1,
        'memory': '4g',
        'walltime': '02:00:00',
        'account': 'EcoGenetics'
    }
    
    spec = f"""
    # Setting conda environments
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate python_phylo
    
    # Writing job information to standard output
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    # Setting working directory
    mkdir -p {work_path}/{spid}
    mkdir -p {log_path}/{spid}
    
    cd {work_path}/{spid}
    
    # Generate FST-based phylogeny
    python {script_path}/fst_phylo_generator.py \\
        {inputs["fst_file"]} \\
        -o {outputs["newick"]} \\
        -m {method}
    
    echo "FINISH: $(date)"
    jobinfo $SLURM_JOBID
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def prepare_iqtree_single_rep(work_path: str, script_path: str, log_path: str, spid: str, 
                             replicate: int, output_name: str, keep_original_names: bool = True):
    """Prepare IQ-TREE inputs for a single bootstrap replicate."""
    
    rep_str = f"{replicate:03d}"
    
    inputs = {
        "fst_tree": f"{work_path}/{spid}/{spid}_fst_phylo.newick",
        "fasta_file": f"{work_path}/{spid}/bootstrap_sequences/{output_name}_rep{rep_str}.fasta",
        "bootstrap_log": f"{log_path}/{spid}/{spid}_bootstrap_sampling.DONE",
        "fst_log": f"{log_path}/{spid}/{spid}_fst_phylo.DONE"
    }
    
    outputs = {
        "topology": f"{work_path}/{spid}/iqtree_inputs/{spid}_rep{rep_str}_topology.newick",
        "sequences": f"{work_path}/{spid}/iqtree_inputs/{spid}_rep{rep_str}_sequences.fasta",
        "mapping": f"{work_path}/{spid}/iqtree_inputs/{spid}_rep{rep_str}_name_mapping.txt",
        "log": f"{log_path}/{spid}/{spid}_iqtree_prep_rep{rep_str}.DONE"
    }
    
    options = {
        'cores': 1,
        'memory': '4g',
        'walltime': '02:00:00',
        'account': 'EcoGenetics'
    }
    
    keep_names_flag = "--keep-original-names" if keep_original_names else ""
    
    spec = f"""
    # Setting conda environments
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate python_phylo
    
    # Writing job information to standard output
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    echo "Processing replicate {rep_str}"
    
    # Setting working directory
    mkdir -p {work_path}/{spid}/iqtree_inputs
    mkdir -p {log_path}/{spid}
    
    cd {work_path}/{spid}/iqtree_inputs
    
    # Check if input files exist
    if [[ ! -f "{inputs["fasta_file"]}" ]]; then
        echo "Error: FASTA file not found: {inputs["fasta_file"]}"
        exit 1
    fi
    
    if [[ ! -f "{inputs["fst_tree"]}" ]]; then
        echo "Error: FST tree file not found: {inputs["fst_tree"]}"
        exit 1
    fi
    
    # Run prepare_phylo_inputs for this replicate
    python {script_path}/prepare_phylo_inputs.py \\
        "{inputs["fst_tree"]}" \\
        "{inputs["fasta_file"]}" \\
        -p "{spid}_rep{rep_str}" \\
        {keep_names_flag}
    
    echo "FINISH: $(date)"
    jobinfo $SLURM_JOBID
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def run_iqtree3_single_rep(work_path: str, script_path: str, log_path: str, spid: str, 
                          replicate: int, model: str = "MFP+G", threads: int = 1):
    """Run IQ-TREE3 on a single bootstrap replicate with constraint tree."""
    
    rep_str = f"{replicate:03d}"
    
    inputs = {
        "topology": f"{work_path}/{spid}/iqtree_inputs/{spid}_rep{rep_str}_topology.newick",
        "sequences": f"{work_path}/{spid}/iqtree_inputs/{spid}_rep{rep_str}_sequences.fasta",
        "prep_log": f"{log_path}/{spid}/{spid}_iqtree_prep_rep{rep_str}.DONE"
    }
    
    outputs = {
        "treefile": f"{work_path}/{spid}/iqtree_results/{spid}_rep{rep_str}.treefile",
        "log": f"{log_path}/{spid}/{spid}_iqtree_rep{rep_str}.DONE"
    }
    
    options = {
        'cores': threads,
        'memory': '24g',
        'walltime': '12:00:00',
        'account': 'EcoGenetics'
    }
    
    spec = f"""
    # Setting conda environments
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate python_phylo
    
    # Writing job information to standard output
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    echo "Running IQ-TREE3 on replicate {rep_str}"
    
    # Setting working directory
    mkdir -p {work_path}/{spid}/iqtree_results
    mkdir -p {log_path}/{spid}
    
    cd {work_path}/{spid}/iqtree_results
    
    # Check if input files exist
    if [[ ! -f "{inputs["sequences"]}" ]]; then
        echo "Error: FASTA file not found: {inputs["sequences"]}"
        exit 1
    fi
    
    if [[ ! -f "{inputs["topology"]}" ]]; then
        echo "Error: Topology file not found: {inputs["topology"]}"
        exit 1
    fi
    
    # Run IQ-TREE3
    echo "Running: iqtree3 -s {inputs["sequences"]} -te {inputs["topology"]} -m {model} -pre {spid}_rep{rep_str} -nt {threads}"
    
    iqtree3 \\
        -s "{inputs["sequences"]}" \\
        -te "{inputs["topology"]}" \\
        -m {model} \\
        -pre {spid}_rep{rep_str} \\
        -nt {threads}
    
    # Check if output was created
    if [[ ! -f "{spid}_rep{rep_str}.treefile" ]]; then
        echo "Error: IQ-TREE3 did not produce expected output file"
        exit 1
    fi
    
    echo "FINISH: $(date)"
    jobinfo $SLURM_JOBID
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
