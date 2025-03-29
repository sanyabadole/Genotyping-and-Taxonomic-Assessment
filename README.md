# Bacterial Genome Analysis Pipeline

A comprehensive pipeline for bacterial genome analysis including taxonomic classification, strain typing, quality assessment, and contamination detection.

## Overview

This pipeline integrates multiple bioinformatics tools to provide complete analysis of bacterial genomes, from taxonomic classification to quality assessment and contamination detection. The workflow automates the process of analyzing multiple genome assemblies, generating comprehensive reports in various formats.

## Team Members

This pipeline was developed by:
- Sanya Badole (sbadole6@gatech.edu)
- Ari J Jimenez
- Yaoheng Li
- Kristupas Paulius
- Rajan Ranjeet Tidke

## Features

- **Genus-level and Species level classification**: Combines MLST and FastANI analysis
- **Quality metrics**: CheckM for completeness/contamination assessment
- **Contamination detection**: Kraken2-based contig-by-contig taxonomy
- **Automated reporting**: HTML/TSV/CSV summaries of all results

## Requirements

- **CPU**: x86_64 architecture recommended
- **OS**: Linux (Tested on Ubuntu 22.04 LTS)
- **Memory**: Minimum 16GB RAM (32GB recommended for large datasets)
- **Storage**: 50GB+ free space for temporary files

## Installation

### Clone the Repository

```bash
git clone https://github.com/sanyabadole/bacterial-analysis-pipeline.git
cd bacterial-analysis-pipeline
```

### Install Miniconda

If you don't have Miniconda already installed:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda init bash
exec bash
```

### Tools Installation

Create and activate the conda environment with all required dependencies:

```bash
# Create environment from the provided environment.yml file
conda env create -f environment.yml

# Activate the environment
conda activate bacterial_analysis_env
```

The environment will automatically install the following dependencies:
- mlst >= 3.0
- checkm-genome >= 1.2
- fastani >= 1.3
- kraken2 >= 2.1
- taxonkit >= 0.15

### Database Installation

#### 1. Kraken2 Database

```bash
# Create directory for databases
mkdir -p databases/kraken2

# Download standard Kraken2 database
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20241228.tar.gz -P databases/kraken2/
tar -xzvf databases/kraken2/k2_standard_16gb_20241228.tar.gz -C databases/kraken2/
```

#### 2. CheckM Database

```bash
# Download and set up CheckM database
mkdir -p databases/checkm
cd databases/checkm
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xzvf checkm_data_2015_01_16.tar.gz
checkm data setRoot $(pwd)
cd ../../
```

## Directory Structure

Create the following directory structure for your project:

```
bacterial-analysis-pipeline/
├── databases/                  # Contains all reference databases
│   ├── kraken2/
│   └── checkm/
├── scripts/                    # Pipeline scripts
│   └── pipeline.sh            # Main pipeline script
├── test_data/                  # Test datasets
│   ├── assemblies/
│   └── gffs/
├── environment.yml             # Conda environment specification
└── README.md                   # This file
```

## Usage

```bash
./scripts/pipeline.sh \
  /path/to/assemblies \  # Directory with .fasta files
  /path/to/gffs \        # Directory with .gff files
  /path/to/output \      # Output directory
  16 \                   # Number of CPU threads
  /path/to/kraken_db     # Kraken2 database directory
```

### Example Commands

Small dataset test:
```bash
./scripts/pipeline.sh test_data/assemblies test_data/gffs test_results 4 databases/kraken2/standard
```

Large dataset analysis:
```bash
./scripts/pipeline.sh \
  /mnt/large_data/assemblies \
  /mnt/large_data/gffs \
  /mnt/large_data/results \
  32 \
  databases/kraken2/standard
```

## Output Structure

```
results/
├── final_summary/             # Comprehensive reports
│   ├── master_summary.tsv
│   ├── report.html
│   └── master_summary.csv
├── genotyping/
│   ├── genus_level/           
│   └── species_level/         # MLST & FastANI results
├── quality_assessment/        # CheckM outputs
└── taxonomy/                  # Kraken2 contamination reports
```

## Pipeline Steps

1. **Genus Classification**
   - Input: .gff annotation files
   - Output: genotyping/genus_level/*.txt

2. **Species Identification**
   - Tools: MLST + FastANI
   - Output: genotyping/species_level/combined_results.tsv

3. **Quality Control**
   - Tool: CheckM
   - Output: quality_assessment/checkm_quality.tsv

4. **Contamination Check**
   - Tool: Kraken2 (requires database)
   - Output: taxonomy/summaries/*_contamination_report.txt

5. **Final Report**
   - Combined HTML/TSV/CSV outputs

## Troubleshooting

### Common Issues

1. **Insufficient Memory**
   - Error: Process killed due to memory allocation failure
   - Solution: Increase available RAM or reduce thread count

2. **Database Path Errors**
   - Error: Cannot find Kraken2 database
   - Solution: Ensure database paths are absolute and correctly specified

3. **Missing Dependencies**
   - Error: Command not found
   - Solution: Ensure conda environment is activated with `conda activate bacterial-analysis`

4. **Time Usage Error**
   - Error: illegal option --f
   - Solution: This usually happens on macOS system with no gnu-time installed, simply run `brew install gnu-time` and restart the pipeline. 
