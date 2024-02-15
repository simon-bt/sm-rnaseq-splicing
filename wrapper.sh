#!/usr/bin/env bash
#SBATCH --time=48:00:00
#SBATCH --qos=xlong
#SBATCH --nodes=1
#SBATCH -J sm_rnaseq_splicing
#SBATCH -o sm_rnaseq_splicing_%A.out
#SBATCH --mem-per-cpu=5G

## Load required env modules
module load Mamba
module load snakemake/7.22.0-foss-2022a
module load Singularity/3.10.0

## Remove logs from a previous run
snakemake --snakefile Snakefile clean -c 1

## Execute entire pipeline
snakemake --snakefile Snakefile --slurm --profile profile/ -c 2
# or: snakemake --snakefile Snakefile --slurm --rerun-incomplete --profile profile/ -c 2

## Clean up STAR directory
snakemake --snakefile Snakefile star_cleanup --slurm --profile profile/ -c 1

## Re-run alternative splicing and sashimi plots with modified config (groupA and groupB)
# snakemake --snakefile Snakefile vast_compare visualize_data prepare_sashimis plot_sashimis \
#   --slurm --profile profile/ -c 2
