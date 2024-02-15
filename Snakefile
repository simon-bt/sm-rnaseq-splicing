# Title: sm-rnaseq-splicing
# Description: Pipeline for alternative splicing analysis with vast-tools
# Author: S. Bajew
# Date: February 2024

#TODO
# 1) Add report generation
# 2) Extend process module with a function for pairwise comparison between two groups
# 3) Add a rule for pairwise comparison
# 4) Add figures to visualize pairwise comparisons

from snakemake.utils import min_version
import pandas


# ==========
# Workflow configuration - DO NOT MODIFY
# ==========

min_version('7.22')
configfile: "config/config.yaml"

# Variables for data processing
SAMPLES, = glob_wildcards("data/{sample}_1.fastq.gz")
N_SAMPLES = len(SAMPLES)
PAIRED_END = config['paired_end']
READ_LENGTH = config['read_length']
COMBINE_OPTS = config['combine_opts']
SPECIES = config['species']
GROUP_A = config['groupA']
GROUP_B = config['groupB']
PAIRED = config['paired']
MIN_DPSI = config['min_dpsi']
MIN_RANGE = config['min_range']
CONFIG = pandas.read_csv(config['samples_path'],sep="\t")
SAMPLES_A = ','.join(CONFIG[CONFIG['groupName'] == GROUP_A]['sampleName'].tolist())
SAMPLES_B = ','.join(CONFIG[CONFIG['groupName'] == GROUP_B]['sampleName'].tolist())
SAMPLES_LIST = CONFIG[CONFIG['groupName'] == GROUP_A]['sampleName'].tolist() + \
               CONFIG[CONFIG['groupName'] == GROUP_B]['sampleName'].tolist()

# Variables for data visualization
COLOR_GROUPA = config['color_groupA']
COLOR_GROUPB = config['color_groupB']
COLOR_REG = config['color_reg']
COLOR_UNREG = config['color_unreg']

# Variables for sashimi plot generation
LOWER_LIMIT = config['lower_limit']
UPPER_LIMIT = config['upper_limit']
SASHIMI_CONTROL = config['sashimi_control']
SASHIMI_CONTRAST = config['sashimi_contrast']
SASHIMI_MINCOV = config['sashimi_min_cov']
SASHIMI_AGGFUNC = config['sashimi_agg_func']


def input_fastq(wildcards):
    """
    Helper function that matches paired-end FASTQ files.

    :param wildcards: Access Snakemake wildcards
    :return: Dictionary of FASTQ files per sample.
    """
    if PAIRED_END:
        return dict(
            zip(['r1', 'r2'],
                expand("data/{sample}_{group}.fastq.gz",group=['1', '2'],**wildcards))
        )
    else:
        return {'r1': "data/{sample}_1.fastq.gz".format(**wildcards)}


# ==========
# Main
# ==========

rule all:
    input:
        f"results/vast_out/{GROUP_B}-vs-{GROUP_A}_dPSI-{MIN_DPSI}.tab",
        f"results/{GROUP_B}_vs_{GROUP_A}/data/DataDiff-{GROUP_B}-vs-{GROUP_A}_dPSI-{MIN_DPSI}.tsv",
        f"results/{GROUP_B}_vs_{GROUP_A}/figures/sashimis/done.txt"


# ==========
# Fetch singularity image
# ==========

rule fetch_singularity:
    output:
        "singularities/ggsashimi.sif"
    envmodules:
        "Singularity/3.10.0"
    threads: 1
    priority: 12
    message: "Fetching singularity image for ggsashimi ..."
    script:
        "scripts/fetch_singularity.sh"


# ==========
# Download files
# ==========

rule download_files:
    output:
        "files/genome.fa",
        "files/annotation.gtf",
        "files/prot_impact.tab",
        "files/event_mapping.tab",
        "files/event_info.tab"
    threads: 4,
    priority: 11
    message: "Downloading required files ..."
    script:
        "scripts/download_files.sh"


# ==========
# Align with vast-tools
# ==========

rule vast_align:
    input:
        unpack(input_fastq)
    output:
        "results/vast_out/to_combine/{sample}.IR.summary_v2.txt",
        "results/vast_out/to_combine/{sample}.IR2",
        "results/vast_out/to_combine/{sample}.MULTI3X",
        "results/vast_out/to_combine/{sample}.eej2",
        "results/vast_out/to_combine/{sample}.exskX",
        "results/vast_out/to_combine/{sample}.info",
        "results/vast_out/to_combine/{sample}.micX"
    envmodules:
        "VAST-TOOLS/2.5.1"
    params:
        paired_end=PAIRED_END,
        species=SPECIES
    threads: 6
    priority: 10
    message: "Aligning {wildcards.sample} RNA-seq reads with vast-tools ..."
    log:
        stderr="logs/vt-align-{sample}.log"
    script:
        "scripts/vast_align.sh"


# ==========
# Combine vast-tools output
# ==========

rule vast_combine:
    input:
        expand('results/vast_out/to_combine/{sample}.IR.summary_v2.txt',sample=SAMPLES),
        expand('results/vast_out/to_combine/{sample}.IR2',sample=SAMPLES),
        expand('results/vast_out/to_combine/{sample}.MULTI3X',sample=SAMPLES),
        expand('results/vast_out/to_combine/{sample}.eej2',sample=SAMPLES),
        expand('results/vast_out/to_combine/{sample}.exskX',sample=SAMPLES),
        expand('results/vast_out/to_combine/{sample}.info',sample=SAMPLES),
        expand('results/vast_out/to_combine/{sample}.micX',sample=SAMPLES)
    output:
        f"results/vast_out/INCLUSION_LEVELS_FULL-{SPECIES}-{N_SAMPLES}.tab"
    envmodules:
        "VAST-TOOLS/2.5.1"
    params:
        opts=COMBINE_OPTS,
        species=SPECIES
    threads: 5
    priority: 9
    message: "Combining vast-tools output into inclusion table ..."
    log:
        stderr="logs/vt-combine.log"
    script:
        "scripts/vast_combine.sh"


# ==========
# Generate alternative splicing data
# ==========

rule vast_compare:
    input:
        f"results/vast_out/INCLUSION_LEVELS_FULL-{SPECIES}-{N_SAMPLES}.tab"
    output:
        f"results/vast_out/AllEvents-{GROUP_B}-vs-{GROUP_A}_dPSI-{MIN_DPSI}.tab",
        f"results/vast_out/{GROUP_B}-vs-{GROUP_A}_dPSI-{MIN_DPSI}.tab"
    envmodules:
        "VAST-TOOLS/2.5.1"
    params:
        species=SPECIES,
        min_dpsi=MIN_DPSI,
        min_range=MIN_RANGE,
        paired=PAIRED,
        groupA=GROUP_A,
        groupB=GROUP_B,
        samplesA=SAMPLES_A,
        samplesB=SAMPLES_B
    threads: 1
    priority: 8
    message: "Generating alternative splicing data ..."
    script:
        "scripts/vast_compare.sh"


# ==========
# Process and visualize alternative splicing data
# ==========

rule process_data:
    input:
        f"results/vast_out/INCLUSION_LEVELS_FULL-{SPECIES}-{N_SAMPLES}.tab",
        f"results/vast_out/AllEvents-{GROUP_B}-vs-{GROUP_A}_dPSI-{MIN_DPSI}.tab",
        f"results/vast_out/{GROUP_B}-vs-{GROUP_A}_dPSI-{MIN_DPSI}.tab",
        "files/prot_impact.tab",
        "files/event_mapping.tab"
    output:
        f"results/{GROUP_B}_vs_{GROUP_A}/data/DataDiff-{GROUP_B}-vs-{GROUP_A}_dPSI-{MIN_DPSI}.tsv",
        touch(f"results/{GROUP_B}_vs_{GROUP_A}/figures/splicing/done.txt"),
        touch(f"results/{GROUP_B}_vs_{GROUP_A}/figures/regulated_events/done.txt")
    conda:
        "envs/config.yaml"
    params:
        samples_a=SAMPLES_A,
        samples_b=SAMPLES_B,
        groupA=GROUP_A,
        groupB=GROUP_B,
        min_dpsi=MIN_DPSI,
        samples=SAMPLES_LIST,
        color_groupA=COLOR_GROUPA,
        color_groupB=COLOR_GROUPB,
        color_reg=COLOR_REG,
        color_unreg=COLOR_UNREG
    threads: 1
    priority: 7
    message: "Processing and visualizing alternative splicing data ..."
    script:
        "scripts/process_data.py"


# ==========
# Generate STAR index
# ==========

rule star_index:
    input:
        "files/genome.fa",
        "files/annotation.gtf"
    output:
        directory("results/star_out/star_index")
    conda:
        "envs/config.yaml"
    params:
        read_length=READ_LENGTH
    threads: 5
    priority: 6
    message: "Generating STAR index ..."
    script:
        "scripts/star_index.sh"


# ==========
# Map with STAR
# ==========

rule star_bams:
    input:
        unpack(input_fastq),
        "results/star_out/star_index"
    output:
        'results/star_out/{sample}.Aligned.sortedByCoord.out.bam',
        'results/star_out/{sample}.done.txt'
    conda:
        "envs/config.yaml"
    params:
        paired_end=PAIRED_END
    threads: 5
    priority: 5
    message: "Mapping {wildcards.sample} RNA-seq reads with STAR..."
    script:
        "scripts/star_bams.sh"


# ==========
# Prepare configs and data for ggsashimi
# ==========

rule prepare_sashimis:
    input:
        "config/samples.tsv",
        f"results/{GROUP_B}_vs_{GROUP_A}/data/DataDiff-{GROUP_B}-vs-{GROUP_A}_dPSI-{MIN_DPSI}.tsv",
        "files/event_info.tab",
        expand("results/star_out/{sample}.done.txt", sample=SAMPLES)
    output:
        directory(f"results/{GROUP_B}_vs_{GROUP_A}/figures/sashimis"),
        f"results/{GROUP_B}_vs_{GROUP_A}/figures/sashimis/config_bams.tsv",
        f"results/{GROUP_B}_vs_{GROUP_A}/figures/sashimis/config.tab",
        f"results/{GROUP_B}_vs_{GROUP_A}/figures/sashimis/palette.txt"
    params:
        lower_limit=LOWER_LIMIT,
        upper_limit=UPPER_LIMIT,
        sashimi_control=SASHIMI_CONTROL,
        sashimi_contrast=SASHIMI_CONTRAST
    threads: 1
    priority: 4
    message: f"Preparing configuration files for sashimi plots for {GROUP_B}_vs_{GROUP_A} comparison ..."
    script:
        "scripts/prepare_sashimis.py"


# ==========
# Plot sashimis
# ==========

rule plot_sashimis:
    input:
        f"results/{GROUP_B}_vs_{GROUP_A}/figures/sashimis/config_bams.tsv",
        f"results/{GROUP_B}_vs_{GROUP_A}/figures/sashimis/config.tab",
        f"results/{GROUP_B}_vs_{GROUP_A}/figures/sashimis/palette.txt",
        "files/annotation.gtf",
        f"results/{GROUP_B}_vs_{GROUP_A}/figures/sashimis/",
        expand("results/star_out/{sample}.done.txt", sample=SAMPLES)
    output:
        touch(f"results/{GROUP_B}_vs_{GROUP_A}/figures/sashimis/done.txt")
    singularity:
        "singularities/ggsashimi.sif"
    params:
        sashimi_min_cov=SASHIMI_MINCOV,
        sashimi_agg_func=SASHIMI_AGGFUNC
    threads: 2
    priority: 3
    message: f"Generating sashimi plots for {GROUP_B}_vs_{GROUP_A} comparison ..."
    script:
        "scripts/plot_sashimis.sh"


# ==========
# Clean /results/star_out directory
# ==========

rule star_cleanup:
    input:
        "results/star_out"
    shell:
        """
        rm -rf {input}/*_STARpass1 {input}/*._STARgenome \
        {input}/*.SJ.out.tab {input}/*.Log.progress.out \
        {input}/*.Log.out _STARtmp/ Log.out
        """


# ==========
# Clean .snakemake directory
# ==========
rule clean:
    shell: "rm -rf .snakemake/"
