# sm-rnaseq-splicing

Snakemake pipeline for alternative splicing analysis with VAST-TOOLS.

Table of Contents:
- [Getting started](#getting-started)
  - [Requirements](#requirements)
  - [Installing](#installing)
  - [Running sm-rnaseq-splicing](#running-sm-rnaseq-splicing)
- [Output](#output)
- [Authors](#authors)
- [License](#license)

## Getting started

### Requirements

**sm-rnaseq-splicing** requires the following:

- [VAST-TOOLS](https://github.com/vastgroup/vast-tools) and its dependencies (environment module)
- [Snakemake](https://snakemake.readthedocs.io/) >= 7.22 (environment module)
- [Mamba](https://mamba.readthedocs.io) >= 23.1.0 (environment module)
- [Singularity](https://sylabs.io/singularity/) (environment module)
- Python >= 3.6


### Installing

Clone this repo `git clone https://github.com/simon-bt/sm-rnaseq-splicing.git`. 
The repo contains the structure required for the pipeline execution.


### Running sm-rnaseq-splicing

1. Modify **config/samples.tsv** file to provide sample and group information.
2. Modify **config/config.yaml** file to specify the type of reads and parameters used for alternative splicing analysis and downstream data visualisation.
3. Modify **profile/config.yaml** file to better suit your SLURM scheduler configuration and assign per-rule job requirements. 
4. Provide soft-link paths to FASTQ files in **{1|2}_fastq.gz** format in **data/** folder.
5. Modify **wrapper.sh** file to better suit your SLURM scheduler configuration.
6. Modify path  of 'singularity-args' in **profile/config.yaml** to bind root directory using Singularity.

It is recommended to execute a dry run from the top directory _before_ submitting **wrapper.sh** to the scheduler: 
`snakemake --snakefile Snakefile --slurm --profile profile/ --dry-run`

## Output

sm-rnaseq-splicing pipeline uses vast-tools to quantify alternative splicing changes from RNA-sequencing data.
The pipeline aligns sequencing reads to a database for a given species, combines vast-tools output into inclusion table 
and allows to compare two conditions at a time of runtime. The comparisons can be performed for any given two conditions 
upon their specification in the **config/config.yaml** file. Pipeline can be also re-run with a different dPSI threshold.

**vast_compare** and **visualize_data** rules generate:

1) In **results/{GroupB}\_vs\_{GroupA}/data**:

* DataAll-*.tsv files - inclusion levels for splicing events (all and separately) passing coverage and other default criteria. Events passing criteria for the regulation are indicated.
* DataDiff-*.tsv files - more extended inclusion level data for the regulated splicing events (all and separately)

**visualize_data** rule generates figures for microexons (exons shorter than 28 nt), exons (longer than 27 nt) and introns:

2) In **results/{GroupB}\_vs\_{GroupA}/figures/splicing/** (pdf format): 

* Scatter_* - scatter plot of splicing quantification levels between two conditions for a given event type
* DistDPSI - box distribution of dPSI values for the regulated events
* LevelsDiff - violin distribution of splicing quantification levels for the regulated events
* RidgeDPSI - ridge density plot for changes in splicing quantification levels for the regulated events
* ProteinImpact - stacked % bar plot of the protein impact predictions for the regulated events
* HMapLevels_Standardized - heatmap of splicing quantification levels for the regulated events, standardized to 1
* HMapLevels_ZScores - heatmap of splicing quantification levels for the regulated events, in z-score of quantification values
* PercentageDPSI - stacked % bar plot of the proportion of truly regulated events that pass vast-tools compare criteria to all events 
  that pass the coverage filter

3) In **results/{GroupB}\_vs\_{GroupA}/figures/regulated_events/** (pdf format) - per-sample quantification level plot
  for the regulated microexons (exons shorter than 28 nt), exons (longer than 27 nt) and introns

**prepare_sashimis** and **plot_sashimis** rules generate sashimi plots for the regulated microexons (exons shorter than 28 nt), 
exons (longer than 27 nt) and introns:

4) In **results/{GroupB}\_vs\_{GroupA}/figures/sashimis/** (pdf format) - configuration files and sashimi plots 


## Authors

Simon Bajew, PhD (IIS Biodonostia)

## License

This project is under MIT License.