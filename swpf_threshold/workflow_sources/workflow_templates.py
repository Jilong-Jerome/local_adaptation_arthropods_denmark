from gwf import AnonymousTarget
import os

def lr_distribution_analysis_template(
    work_path: str,
    script_path: str,
    log_path: str,
    spid: str,
    input_file: str,
    output_dir: str
):
    """
    Template for analyzing SweepFinder2 likelihood ratio distribution per population.
    
    Args:
        work_path: Working directory path
        script_path: Path to scripts directory
        log_path: Path to log directory
        spid: Species ID
        input_file: Path to SweepFinder2 output file
        output_dir: Output directory for distribution analysis results
    """
    
    inputs = [input_file]
    outputs = [
        f"{output_dir}/LR_summary_statistics.tsv",
        f"{log_path}/{spid}/{spid}_lr_distribution.DONE"
    ]
    
    options = {
        'cores': 1,
        'memory': '8g',
        'walltime': '02:00:00',
        'account': 'EcoGenetics'
    }
    
    spec = f'''
# SweepFinder2 Likelihood Ratio Distribution Analysis for {spid}
echo "Starting LR distribution analysis for {spid}..."
echo "Input file: {input_file}"
echo "Output directory: {output_dir}"

# Create output and log directories
mkdir -p {output_dir}
mkdir -p {log_path}/{spid}

# Activate conda environment
source /home/jilong/miniforge3/etc/profile.d/conda.sh
conda activate python_phylo

# Run the distribution analysis script
cd {work_path}
python3 {script_path}/summarize_lr_distribution.py {input_file} {output_dir}

# Check if analysis completed successfully
if [ $? -eq 0 ]; then
    echo "LR distribution analysis completed successfully for {spid}"
    touch {log_path}/{spid}/{spid}_lr_distribution.DONE
else
    echo "LR distribution analysis failed for {spid}"
    exit 1
fi

# Deactivate conda environment
conda deactivate
'''
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
