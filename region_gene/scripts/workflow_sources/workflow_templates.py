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


def blast_gene_annotation(work_path: str, script_path: str, log_path: str, spid: str,
                         blast_db: str, region_id: str = None, evalue: float = 1e-5):
    """Run BLAST annotation on extracted gene sequences."""
    
    if region_id:
        # Process specific region
        inputs = {
            "fasta_file": f"{work_path}/{spid}/region_genes/{region_id}_genes.fasta",
            "extract_log": f"{log_path}/{spid}/{spid}_extract_genes_all.DONE"
        }
        outputs = {
            "annotation": f"{work_path}/{spid}/gene_annotations/{region_id}_gene_annotations.tsv",
            "log": f"{log_path}/{spid}/{spid}_diamond_{region_id}.DONE"
        }
        region_arg = f"--region-id {region_id}"
    else:
        # Process all regions - this would need to be implemented differently
        # For now, focus on single region processing
        inputs = {
            "gene_dir": f"{work_path}/{spid}/region_genes/",
            "extract_log": f"{log_path}/{spid}/{spid}_extract_genes_all.DONE"
        }
        outputs = {
            "annotation_dir": f"{work_path}/{spid}/gene_annotations/",
            "log": f"{log_path}/{spid}/{spid}_diamond_all.DONE"
        }
        region_arg = ""
    
    options = {
        'cores': 4,
        'memory': '16g',
        'walltime': '12:00:00',
        'account': 'EcoGenetics'
    }
    
    if region_id:
        # Single region processing
        spec = f"""
        # Setting conda environments
        CONDA_BASE=$(conda info --base)
        source $CONDA_BASE/etc/profile.d/conda.sh
        conda activate python_phylo
        
        # Writing job information to standard output
        echo "START: $(date)"
        echo "JobID: $SLURM_JOBID"
        
        # Setting working directory
        mkdir -p {work_path}/{spid}/gene_annotations
        mkdir -p {log_path}/{spid}
        
        cd {work_path}/{spid}/gene_annotations
        
        # Check if input FASTA file exists
        if [[ ! -f "{inputs["fasta_file"]}" ]]; then
            echo "Error: FASTA file not found: {inputs["fasta_file"]}"
            exit 1
        fi
        
        # Run Diamond annotation
        python {script_path}/blast_gene_annotation.py \\
            --fasta-file {inputs["fasta_file"]} \\
            --diamond-db {blast_db}/reference_proteomes.dmnd \\
            --output-file {outputs["annotation"]} \\
            --log-file {work_path}/{spid}/gene_annotations/diamond_{region_id}.log \\
            --evalue {evalue}
        
        echo "FINISH: $(date)"
        jobinfo $SLURM_JOBID
        echo done > {outputs["log"]}
        """
    else:
        # All regions processing
        spec = f"""
        # Setting conda environments
        CONDA_BASE=$(conda info --base)
        source $CONDA_BASE/etc/profile.d/conda.sh
        conda activate python_phylo
        
        # Writing job information to standard output
        echo "START: $(date)"
        echo "JobID: $SLURM_JOBID"
        
        # Setting working directory
        mkdir -p {work_path}/{spid}/gene_annotations
        mkdir -p {log_path}/{spid}
        
        cd {work_path}/{spid}/gene_annotations
        
        # Process all FASTA files in region_genes directory
        for fasta_file in {work_path}/{spid}/region_genes/*_genes.fasta; do
            if [[ -f "$fasta_file" ]]; then
                basename=$(basename "$fasta_file" _genes.fasta)
                echo "Processing $basename..."
                
                python {script_path}/blast_gene_annotation.py \\
                    --fasta-file "$fasta_file" \\
                    --diamond-db {blast_db}/reference_proteomes.dmnd \\
                    --output-file "{work_path}/{spid}/gene_annotations/${{basename}}_gene_annotations.tsv" \\
                    --log-file {work_path}/{spid}/gene_annotations/diamond_${{basename}}.log \\
                    --evalue {evalue}
            fi
        done
        
        echo "FINISH: $(date)"
        jobinfo $SLURM_JOBID
        echo done > {outputs["log"]}
        """
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def extract_all_genes(work_path: str, script_path: str, log_path: str, spid: str,
                     genome_file: str, gff_file: str):
    """Extract all gene sequences from genome using GFF annotation into a single FASTA file."""
    inputs = {
        "genome_file": genome_file, 
        "gff_file": gff_file
    }
    
    outputs = {
        "fasta_file": f"{work_path}/{spid}/all_genes/all_genes.fasta",
        "log": f"{log_path}/{spid}/{spid}_extract_all_genes.DONE"
    }
    
    options = {
        'cores': 2,
        'memory': '16g',
        'walltime': '04:00:00',
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
    mkdir -p {work_path}/{spid}/all_genes
    mkdir -p {log_path}/{spid}
    
    cd {work_path}/{spid}/all_genes
    
    # Check if input files exist
    if [[ ! -f "{inputs["genome_file"]}" ]]; then
        echo "Error: Genome file not found: {inputs["genome_file"]}"
        exit 1
    fi
    
    if [[ ! -f "{inputs["gff_file"]}" ]]; then
        echo "Error: GFF file not found: {inputs["gff_file"]}"
        exit 1
    fi
    
    # Run whole-genome gene extraction
    python {script_path}/extract_all_genes.py \\
        --genome-file {inputs["genome_file"]} \\
        --gff-file {inputs["gff_file"]} \\
        --output-dir {work_path}/{spid}/all_genes \\
        --log-file {work_path}/{spid}/all_genes/gene_extraction.log \\
        --output-filename all_genes.fasta
    
    echo "FINISH: $(date)"
    jobinfo $SLURM_JOBID
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def blast_all_genes_annotation(work_path: str, script_path: str, log_path: str, spid: str,
                               blast_db: str, evalue: float = 1e-5):
    """Run BLAST annotation on all genes."""
    
    inputs = {
        "fasta_file": f"{work_path}/{spid}/all_genes/all_genes.fasta",
        "extract_log": f"{log_path}/{spid}/{spid}_extract_all_genes.DONE"
    }
    outputs = {
        "annotation": f"{work_path}/{spid}/all_gene_annotations/all_genes_annotations.tsv",
        "log": f"{log_path}/{spid}/{spid}_diamond_all_genes.DONE"
    }
    
    options = {
        'cores': 8,
        'memory': '32g',
        'walltime': '24:00:00',
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
    mkdir -p {work_path}/{spid}/all_gene_annotations
    mkdir -p {log_path}/{spid}
    
    cd {work_path}/{spid}/all_gene_annotations
    
    # Check if input FASTA file exists
    if [[ ! -f "{inputs["fasta_file"]}" ]]; then
        echo "Error: FASTA file not found: {inputs["fasta_file"]}"
        exit 1
    fi
    
    # Run Diamond annotation on all genes
    python {script_path}/blast_gene_annotation.py \\
        --fasta-file {inputs["fasta_file"]} \\
        --diamond-db {blast_db}/reference_proteomes.dmnd \\
        --output-file {outputs["annotation"]} \\
        --log-file {work_path}/{spid}/all_gene_annotations/diamond_all_genes.log \\
        --evalue {evalue}
    
    echo "FINISH: $(date)"
    jobinfo $SLURM_JOBID
    echo done > {outputs["log"]}
    """
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
