#!/bin/env python3
from gwf import AnonymousTarget
import os, glob


def extract_region_genes(work_path: str, script_path: str, log_path: str, spid: str,
                        regions_file: str, genome_file: str, gff_file: str, min_length: int = 5000,
                        flank_size: int = 10000, region_id: str = None):
    """Extract gene sequences from selection regions."""
    inputs = {
        "regions_file": regions_file,
        "genome_file": genome_file, 
        "gff_file": gff_file
    }
    
    if region_id:
        # Process specific region
        outputs = {
            "fasta": f"{work_path}/{spid}/region_genes/{region_id}_genes.fasta",
            "log": f"{log_path}/{spid}/{spid}_extract_genes_{region_id}.DONE"
        }
        region_arg = f"--region-id {region_id}"
    else:
        # Process all regions
        outputs = {
            "output_dir": f"{work_path}/{spid}/region_genes/",
            "log": f"{log_path}/{spid}/{spid}_extract_genes_all.DONE"
        }
        region_arg = ""
    
    options = {
        'cores': 1,
        'memory': '8g',
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
    mkdir -p {work_path}/{spid}/region_genes
    mkdir -p {log_path}/{spid}
    
    cd {work_path}/{spid}/region_genes
    
    # Check if input files exist
    if [[ ! -f "{inputs["regions_file"]}" ]]; then
        echo "Error: Regions file not found: {inputs["regions_file"]}"
        exit 1
    fi
    
    if [[ ! -f "{inputs["genome_file"]}" ]]; then
        echo "Error: Genome file not found: {inputs["genome_file"]}"
        exit 1
    fi
    
    if [[ ! -f "{inputs["gff_file"]}" ]]; then
        echo "Error: GFF file not found: {inputs["gff_file"]}"
        exit 1
    fi
    
    # Run gene extraction
    python {script_path}/extract_region_genes.py \\
        --regions-file {inputs["regions_file"]} \\
        --genome-file {inputs["genome_file"]} \\
        --gff-file {inputs["gff_file"]} \\
        --output-dir {work_path}/{spid}/region_genes \\
        --log-file {work_path}/{spid}/region_genes/gene_extraction.log \\
        --min-length {min_length} \\
        --flank-size {flank_size} \\
        {region_arg}
    
    echo "FINISH: $(date)"
    jobinfo $SLURM_JOBID
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)