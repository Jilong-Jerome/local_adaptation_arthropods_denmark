# local_adaptation_arthropods_denmark
local adaptation in Danish arthropods
![til](https://github.com/Jilong-Jerome/local_adaptation_arthropods_denmark/blob/main/all_common_selections.gif)

---

## Bootstrap Phylogeny Sampler

A Python script for generating bootstrap sequences for phylogenetic analysis from pooled population VCF files.

### Overview

`bootstrap_phylo_sampler.py` samples sequences for phylogenetic tree construction from pooled sequencing VCF data. The script randomly samples alleles at each SNP position based on the allele frequencies (read depths) for each population, creating concatenated sequences suitable for phylogenetic analysis.

### Features

- **Random sampling**: Samples alleles based on population-specific allele frequencies from read depths
- **Bootstrap replicates**: Generate multiple bootstrap replicates for phylogenetic confidence estimation
- **Flexible filtering**: Filter by specific scaffolds or genomic regions
- **Concatenated sequences**: Creates one sequence per population containing all selected SNPs
- **Reproducible results**: Set random seed for consistent sampling
- **Organized output**: Specify output directory for result files

### Requirements

- Python 3.x
- Input VCF file with pooled population data containing AD (Allelic Depth) fields

### Usage

```bash
./bootstrap_phylo_sampler.py [options] input.vcf
```

#### Basic Examples

**Generate single concatenated sequence per population:**
```bash
./bootstrap_phylo_sampler.py data/EntNic_sweepfinder_input.vcf -o phylo_sequences.fasta
```

**Generate 100 bootstrap replicates:**
```bash
./bootstrap_phylo_sampler.py data/EntNic_sweepfinder_input.vcf --outdir bootstrap_results -o phylo_sequences -r 100 --seed 42
```

**Filter to specific scaffolds:**
```bash
./bootstrap_phylo_sampler.py data/EntNic_sweepfinder_input.vcf -s HiC_scaffold_1 HiC_scaffold_2 -r 100 --outdir scaffold_analysis
```

**Filter to specific genomic regions:**
```bash
./bootstrap_phylo_sampler.py data/EntNic_sweepfinder_input.vcf --regions HiC_scaffold_1:1000000-2000000 HiC_scaffold_2:500000-1500000 -r 100
```

#### Command Line Options

| Option | Description |
|--------|-------------|
| `vcf` | Input VCF file (can be gzipped) |
| `-o, --output` | Output FASTA filename (default: phylo_sequences.fasta) |
| `-r, --replicates` | Number of bootstrap replicates (default: 1) |
| `-s, --scaffolds` | Specific scaffolds to include (e.g., HiC_scaffold_1 HiC_scaffold_2) |
| `--regions` | Specific regions to include (format: scaffold:start-end or scaffold) |
| `--seed` | Random seed for reproducibility |
| `--outdir` | Output directory for FASTA files |

#### Region Format

Regions can be specified in several formats:
- `scaffold` - Include entire scaffold
- `scaffold:position` - Single position
- `scaffold:start-end` - Range of positions

Examples:
- `HiC_scaffold_1` - Entire scaffold 1
- `HiC_scaffold_1:1000000-2000000` - Region from 1Mb to 2Mb on scaffold 1
- `HiC_scaffold_2:500000` - Single position at 500kb on scaffold 2

### Output

#### Single Replicate
When `--replicates 1` (default), creates one FASTA file:
- `phylo_sequences.fasta` (or specified filename)

#### Multiple Replicates
When `--replicates > 1`, creates numbered FASTA files:
- `phylo_sequences_rep001.fasta`
- `phylo_sequences_rep002.fasta`
- ...
- `phylo_sequences_rep100.fasta`

#### FASTA Format
Each file contains one sequence per population:
```
>EntNic_aaRJ-C225_rep001
ATGCATGCATGCATGCATGC...
>EntNic_aeRoe-C36_rep001
ATCCATGCATGGATGCATGC...
>EntNic_BIJ-C30_rep001
ATGGATGCATGCATGCATGC...
```

### Methodology

1. **Allele Frequency Calculation**: For each SNP position and population, calculates allele frequencies from AD (Allelic Depth) fields
2. **Random Sampling**: Samples ref/alt alleles based on their frequencies in each population
3. **Sequence Construction**: Concatenates sampled alleles across all SNPs to create full sequences
4. **Bootstrap Replicates**: Repeats sampling process for specified number of replicates

#### Missing Data Handling
- Sites with no coverage (AD = 0,0) use reference allele (no variants detected)
- Multi-allelic sites are skipped
- Sites without AD fields are skipped

