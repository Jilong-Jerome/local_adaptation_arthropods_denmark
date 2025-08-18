#!/bin/env python3
from gwf import AnonymousTarget
import os, glob


def parse_input_sweepfinder2(work_path: str, script_path: str, log_path: str, spid: str, vcf: str, maf : float):
    inputs = {"vcf":vcf}
    outputs = {
            "log":f"{log_path}/{spid}/{spid}_{maf}_input_parsing_sweepfinder2.DONE"
            }
    options = {
            'cores':1,
            'memory':'1g',
            'walltime':'1:00:00',
            'account':'EcoGenetics'
        }
    spec = f"""
    # Setting conda environments
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate biopython
    # Writting job information to standard output
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    # Setting working directory
    mkdir -p {work_path}/{spid}/parsed_inputs
    mkdir -p {log_path}/{spid}
    cd {work_path}/{spid}/parsed_inputs
    ######Partly adapted scripts from Fatima#############1_input_files_count_freq_clean.sh#######
    grep -v "#" {inputs["vcf"]} > {spid}_input_vcf_headless.vcf
    cut -f1,2,10- {spid}_input_vcf_headless.vcf | \\
    python {script_path}/2_python_coloninfo_combscript.py | \\
    awk -F '\\t' '{{count0=0; count1=0; for(i=3; i<=NF; i++) {{count0+=gsub(/0/, "", $i); count1+=gsub(/1/, "", $i)}}; print $0, count0, count1}}' | \\
    awk -F'[ /]+' '{{zero=$(NF-1); one=$NF;total = zero+one; zerofreq=zero/total; onefreq=one/total; maf=(zerofreq<onefreq)?zerofreq:onefreq; if (maf>{maf}) print $0;}}' > {spid}_maf_{maf}.tsv
    awk 'NR==FNR{{columns[$1,$2]; next}} ($1,$2) in columns {{print $0}}' {spid}_maf_{maf}.tsv {spid}_input_vcf_headless.vcf > {spid}_maf_{maf}_population_filtered.tsv
    cut -f1,2,10- {spid}_maf_{maf}_population_filtered.tsv > {spid}_maf_{maf}_chr_pos_genoinfo.tsv
    python {script_path}/3_colon_count_freq.py {spid}_maf_{maf}_chr_pos_genoinfo.tsv {spid}_maf_{maf}_chr_pos_count1.tsv {spid}_maf_{maf}_chr_pos_freq.tsv
    #a file for each chromosome: esp useful for sweepfinder.
    mkdir -p {work_path}/{spid}/parsed_inputs/split_chrs
    cd {work_path}/{spid}/parsed_inputs/split_chrs
    awk '{{print > $1"_maf_{maf}_chr_pos_count1.txt"}}' {work_path}/{spid}/parsed_inputs/{spid}_maf_{maf}_chr_pos_count1.tsv
    ######################################################################################
    rm {work_path}/{spid}/parsed_inputs/{spid}_maf_{maf}.tsv
    rm {work_path}/{spid}/parsed_inputs/{spid}_input_vcf_headless.vcf 
    rm {work_path}/{spid}/parsed_inputs/{spid}_maf_{maf}_population_filtered.tsv
    rm {work_path}/{spid}/parsed_inputs/{spid}_maf_{maf}_chr_pos_genoinfo.tsv
    echo "FINISH: $(date)"
    jobinfo $SLURM_JOBID
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)


def parse_input_sweepfinder2_jilong(work_path: str, script_path: str, log_path: str, spid: str, vcf: str, maf : float):
    inputs = {"vcf":vcf}
    outputs = {
            "log":f"{log_path}/{spid}/{spid}_{maf}_input_parsing_sweepfinder2_jilong.DONE",
            "pop_list":f"{work_path}/{spid}/parsed_inputs/{spid}_maf{maf}_sweepfinder_pop_header.tsv",
            "chr_list":f"{work_path}/{spid}/parsed_inputs/{spid}_maf{maf}_sweepfinder_chr_list.txt"
            }
    options = {
            'cores':1,
            'memory':'8g',
            'walltime':'12:00:00',
            'account':'EcoGenetics'
        }
    spec = f"""
    # Setting conda environments
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate biopython
    # Writting job information to standard output
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    # Setting working directory
    mkdir -p {work_path}/{spid}/parsed_inputs
    mkdir -p {log_path}/{spid}
    cd {work_path}/{spid}/parsed_inputs
    python {script_path}/parse_vcf_to_sweepfinder2_input.py {inputs["vcf"]} {maf} {spid}
    mkdir -p {work_path}/{spid}/parsed_inputs/split_chrs/jilong
    cd {work_path}/{spid}/parsed_inputs/split_chrs/jilong
    awk '{{print > $1"_maf_{maf}_chr_pos_count1.txt"}}' {work_path}/{spid}/parsed_inputs/{spid}_maf{maf}_sweepfinder_input.tsv
    find *chr_pos_count1.txt| sed 's/_maf_.*//' > {work_path}/{spid}/parsed_inputs/{spid}_maf{maf}_sweepfinder_chr_list.txt
    echo "FINISH: $(date)"
    jobinfo $SLURM_JOBID
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def split_sweefinder_input(work_path: str, script_path: str, log_path: str, spid: str, maf : float, chrom: chr,chrom_id: chr, pop: chr, pop_col: int, pop_size: int):
    inputs = {
            "parse_done":f"{log_path}/{spid}/{spid}_{maf}_input_parsing_sweepfinder2_jilong.DONE",
            }
    outputs = {
            "log":f"{log_path}/{spid}/{spid}_{chrom_id}_{pop}_{maf}_split.DONE",
            "pop_input":f"{work_path}/{spid}/parsed_inputs/splits/{spid}_maf{maf}_{chrom_id}_{pop}.txt"
            }
    options = {
            'cores':1,
            'memory':'1g',
            'walltime':'1:00:00',
            'account':'EcoGenetics'
        }
    spec = f"""
    # Setting conda environments
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate biopython
    # Writting job information to standard output
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    # Setting working directory
    mkdir -p {work_path}/{spid}/parsed_inputs/splits
    mkdir -p {log_path}/{spid}
    cd {work_path}/{spid}/parsed_inputs/splits
    cut -f 2,{pop_col} {work_path}/{spid}/parsed_inputs/split_chrs/jilong/'{chrom}_maf_{maf}_chr_pos_count1.txt' > {chrom_id}_{pop}_maf_{maf}_pos_count1.tsv
    awk '{{print $0 "\\t{pop_size}"}}' {chrom_id}_{pop}_maf_{maf}_pos_count1.tsv |sed "s/\\t/ /g" > {spid}_maf{maf}_{chrom_id}_{pop}.data
    echo "position\tx\tn" > {spid}_{chrom_id}_{pop}.header
    cat {spid}_{chrom_id}_{pop}.header {spid}_maf{maf}_{chrom_id}_{pop}.data > {spid}_maf{maf}_{chrom_id}_{pop}.txt
    rm {chrom_id}_{pop}_maf_{maf}_pos_count1.tsv
    rm {spid}_maf{maf}_{chrom_id}_{pop}.data
    rm {spid}_{chrom_id}_{pop}.header
    echo "FINISH: $(date)"
    jobinfo $SLURM_JOBID
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def sweepfinder_chr_pop_spectrum(work_path: str, script_path: str, log_path: str, spid: str, maf : float, chrom: chr, chrom_id: chr, pop: chr):
    inputs = {
            "split_done":f"{log_path}/{spid}/{spid}_{chrom_id}_{pop}_{maf}_split.DONE"
        }
    outputs = {
            "log":f"{log_path}/{spid}/{spid}_{chrom_id}_{pop}_{maf}_sweepfinder2_spectrum.DONE",
            "specturm":f"{work_path}/{spid}/sweepfinder2_output/splits/{chrom_id}/{pop}/{spid}_{pop}_{chrom_id}_maf{maf}_spectrum.txt"
        }
    options = {
            'cores':1,
            'memory':'4g',
            'walltime':'12:00:00',
            'account':'EcoGenetics'
        }
    spec = f"""
    # Setting conda environments
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate sweepfinder2
    # Writting job information to standard output
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    # Setting working directory
    mkdir -p {work_path}/{spid}/sweepfinder2_output/splits/{chrom_id}/{pop}
    mkdir -p {log_path}/{spid}
    cd {work_path}/{spid}/sweepfinder2_output/splits/{chrom_id}/{pop}
    SweepFinder2 -f {work_path}/{spid}/parsed_inputs/splits/{spid}_maf{maf}_{chrom_id}_{pop}.txt {spid}_{pop}_{chrom_id}_maf{maf}_spectrum.txt
    echo "FINISH: $(date)"
    jobinfo $SLURM_JOBID
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)


def sweepfinder_chr_pop_scan(work_path: str, script_path: str, log_path: str, spid: str, maf : float, chrom: chr, chrom_id: chr, pop: chr, grid_size: int, option: int):
    '''
    option 1: test selection with a grid size (interval between tested site)
    option 2: test all variant sites
    '''
    if option == 1:
        scan_cmd = f"SweepFinder2 -lg {grid_size} {spid}_{chrom_id}_maf{maf}_w{grid_size}_o{option}_pop_input.temp {spid}_{chrom_id}_maf{maf}_w{grid_size}_o{option}_pop_spec.temp {spid}_{pop}_{chrom_id}_maf{maf}_w{grid_size}_o{option}_sweepfinder2_output.txt"
    elif option == 2:
        scan_cmd = f"SweepFinder2 -lu {spid}_{chrom_id}_maf{maf}_w{grid_size}_o{option}_pop_input.grid {spid}_{chrom_id}_maf{maf}_w{grid_size}_o{option}_pop_input.temp {spid}_{chrom_id}_maf{maf}_w{grid_size}_o{option}_pop_spec.temp {spid}_{pop}_{chrom_id}_maf{maf}_w{grid_size}_o{option}_sweepfinder2_output.txt"
    inputs = {
            "split_done":f"{log_path}/{spid}/{spid}_{chrom_id}_{pop}_{maf}_split.DONE",
            "spectrum_done":f"{log_path}/{spid}/{spid}_{chrom_id}_{pop}_{maf}_sweepfinder2_spectrum.DONE",
            "pop_input":f"{work_path}/{spid}/parsed_inputs/splits/{spid}_maf{maf}_{chrom_id}_{pop}.txt",
            "pop_specturm":f"{work_path}/{spid}/sweepfinder2_output/splits/{chrom_id}/{pop}/{spid}_{pop}_{chrom_id}_maf{maf}_spectrum.txt"
        }
    outputs = {
            "log":f"{log_path}/{spid}/{spid}_{chrom_id}_{pop}_{maf}_w{grid_size}_o{option}_sweepfinder2.DONE"
        }
    options = {
            'cores':1,
            'memory':'8g',
            'walltime':'72:00:00',
            'account':'EcoGenetics'
        }
    spec = f"""
    # Setting conda environments
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate sweepfinder2
    # Writting job information to standard output
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    # Setting working directory
    mkdir -p {work_path}/{spid}/sweepfinder2_output/splits/{chrom_id}/{pop}
    mkdir -p {log_path}/{spid}
    cd {work_path}/{spid}/sweepfinder2_output/splits/{chrom_id}/{pop}
    # prepare pop_input
    cp {inputs["pop_input"]} {spid}_{chrom_id}_maf{maf}_w{grid_size}_o{option}_pop_input.temp
    cp {inputs["pop_specturm"]} {spid}_{chrom_id}_maf{maf}_w{grid_size}_o{option}_pop_spec.temp
    cut -f 1 {spid}_{chrom_id}_maf{maf}_w{grid_size}_o{option}_pop_input.temp| sed '1d' |cut -f 1 -d " " > {spid}_{chrom_id}_maf{maf}_w{grid_size}_o{option}_pop_input.grid
    {scan_cmd}
    rm {spid}_{chrom_id}_maf{maf}_w{grid_size}_o{option}_pop_input.temp
    rm {spid}_{chrom_id}_maf{maf}_w{grid_size}_o{option}_pop_spec.temp
    sed '1d' {spid}_{pop}_{chrom_id}_maf{maf}_w{grid_size}_o{option}_sweepfinder2_output.txt| \\
    awk '{{print $0 "\\t{spid}\\t{pop}\\t{chrom_id}"}}' > {spid}_{pop}_{chrom_id}_maf{maf}_w{grid_size}_o{option}_swpf2_noheader.tsv
    echo "FINISH: $(date)"
    jobinfo $SLURM_JOBID
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def sweepfinder_concat_res(work_path: str, script_path: str, log_path: str, spid: str,  maf : float, grid_size: int, option: int, chr_pos_logs: list):
    inputs = chr_pos_logs
    outputs = {"log":f"{log_path}/{spid}/{spid}_{maf}_w{grid_size}_o{option}_sweepfinder2_combined.DONE",
            "result":f"{work_path}/{spid}/sweepfinder2_output/{spid}_maf{maf}_w{grid_size}_o{option}_swpf2.tsv"}
    options = {
            'cores':1,
            'memory':'4g',
            'walltime':'12:00:00',
            'account':'EcoGenetics'
        }
    spec = f"""
    # Setting conda environments
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate sweepfinder2
    # Writting job information to standard output
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID":
    # Setting working directory
    mkdir -p {work_path}/{spid}/sweepfinder2_output
    mkdir -p {log_path}/{spid}
    cd {work_path}/{spid}/sweepfinder2_output
    cat {work_path}/{spid}/sweepfinder2_output/splits/*/*/{spid}_*_maf{maf}_w{grid_size}_o{option}_swpf2_noheader.tsv > {spid}_maf{maf}_w{grid_size}_o{option}_swpf2_noheader.tsv
    (echo "location\tLR\talpha\tsp\tpop\tchrom"; cat {spid}_maf{maf}_w{grid_size}_o{option}_swpf2_noheader.tsv) > {spid}_maf{maf}_w{grid_size}_o{option}_swpf2.tsv
    rm {spid}_maf{maf}_w{grid_size}_o{option}_swpf2_noheader.tsv
    echo "FINISH: $(date)"
    jobinfo $SLURM_JOBID
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def sweepfinder_viz_res(work_path: str, script_path: str, log_path: str, spid: str,  maf : float, grid_size: int, option: int):
    inputs = {
            "concat_done":f"{log_path}/{spid}/{spid}_{maf}_w{grid_size}_o{option}_sweepfinder2_combined.DONE",
            "viz_data":f"{work_path}/{spid}/sweepfinder2_output/{spid}_maf{maf}_w{grid_size}_o{option}_swpf2.tsv"
            }
    outputs = {
            "log":f"{log_path}/{spid}/{spid}_{maf}_w{grid_size}_o{option}_sweepfinder2_viz.DONE"
            }
    options = {
            'cores':1,
            'memory':'24g',
            'walltime':'12:00:00',
            'account':'EcoGenetics'
        }
    spec = f"""
    # Setting conda environments
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate Rtidy
    # Writting job information to standard output
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    # Setting working directory
    mkdir -p {work_path}/{spid}/sweepfinder2_results_viz
    mkdir -p {log_path}/{spid}
    cd {work_path}/{spid}/sweepfinder2_results_viz
    echo "{script_path}/viz_chr_pos_res.R"
    echo "Rscript {script_path}/viz_chr_pos_res.R {inputs["viz_data"]} {spid} {maf} {grid_size} {option}"
    Rscript {script_path}/viz_chr_pos_res.R {inputs["viz_data"]} {spid} {maf} {grid_size} {option}
    mkdir -p {work_path}/{spid}/sweepfinder2_results_viz/sinlge
    mv *_single.pdf {work_path}/{spid}/sweepfinder2_results_viz/sinlge
    mkdir -p {work_path}/{spid}/sweepfinder2_results_viz/per_chrom
    mv *_all_pops.pdf {work_path}/{spid}/sweepfinder2_results_viz/per_chrom
    mkdir -p {work_path}/{spid}/sweepfinder2_results_viz/per_pop
    mv *_all_chroms.pdf {work_path}/{spid}/sweepfinder2_results_viz/per_pop
    echo "FINISH: $(date)"
    jobinfo $SLURM_JOBID
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)


