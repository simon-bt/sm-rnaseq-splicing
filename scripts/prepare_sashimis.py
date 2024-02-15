import pandas
import os

# DO NOT MODIFY - Variables for data preparation
SAMPLES = snakemake.input[0]
DF_DIFF = snakemake.input[1]
COORDS = snakemake.input[2]
OUTDIR = snakemake.output[0]
LOWER_LIMIT = snakemake.params['lower_limit']
UPPER_LIMIT = snakemake.params['upper_limit']
SASHIMI_CONTROL = snakemake.params['sashimi_control']
SASHIMI_CONTRAST = snakemake.params['sashimi_contrast']
DIR = os.path.abspath(
    os.getcwd()
)

# Prepare config file for BAM files
SAMPLES_CONFIG = pandas.read_csv(SAMPLES, sep='\t')
SAMPLES_CONFIG['path'] = SAMPLES_CONFIG['sampleName']. \
    apply(lambda x: f"{DIR}/results/star_out/{x}.Aligned.sortedByCoord.out.bam")
SAMPLES_CONFIG[['sampleName', 'path', 'groupName']]. \
    to_csv(f"{OUTDIR}/config_bams.tsv", sep='\t', index=False, header=None)

# Get EventIDs and reference coordinates for regulated events
DIFF_EVENTS = pandas. \
    read_csv(DF_DIFF, sep="\t"). \
    set_index('EVENT'). \
    filter(regex='EX|INT', axis=0). \
    reset_index()['EVENT'].tolist()
EVENT_INFO = pandas.\
    read_csv(COORDS, sep='\t')[['EVENT', 'REF_CO', 'CO_A']]
VAST_DIFF_EVENTS = EVENT_INFO[
    EVENT_INFO['EVENT'].isin(DIFF_EVENTS)].\
        rename(columns={'CO_A':'COORD'})

# Process coordinates for regulated exons
EXONS = VAST_DIFF_EVENTS.\
    set_index('EVENT').\
    filter(regex='EX', axis=0).\
    reset_index()
EXONS[['CHR', 'UPSTREAM', 'START', 'END', 'DOWNSTREAM', 'STRAND']] = EXONS['REF_CO'].\
    str.split(':|,|-', expand=True, n=5)
EXONS_pos = EXONS.query('STRAND == \'+\'')
EXONS_neg = EXONS.query('STRAND == \'-\'').\
    rename(columns={'UPSTREAM': 'DOWNSTREAM', 'DOWNSTREAM': 'UPSTREAM'})

# Process coordinates for regulated introns
INTRONS = VAST_DIFF_EVENTS.\
    set_index('EVENT').\
    filter(regex='INT', axis=0).\
    reset_index()
INTRONS[['CHR', 'UPSTREAM', 'START', 'END', 'DOWNSTREAM', 'STRAND']] = INTRONS['REF_CO'].\
    str.split(':|,|=|-', expand=True, n=5)
INTRONS_pos = INTRONS.query('STRAND == \'+\'')
INTRONS_neg = INTRONS.query('STRAND == \'-\'').\
    rename(columns={'UPSTREAM': 'DOWNSTREAM', 'DOWNSTREAM': 'UPSTREAM'})

# Prepare coordinates with new bounds
CONFIG = pandas.concat([
    EXONS_pos,
    EXONS_neg,
    INTRONS_pos,
    INTRONS_neg]
)
CONFIG['UPSTREAM'] = CONFIG['UPSTREAM'].astype('int')
CONFIG['DOWNSTREAM'] = CONFIG['DOWNSTREAM'].astype('int')
CONFIG.loc[:, 'UPSTREAM'] = CONFIG['UPSTREAM'] - LOWER_LIMIT
CONFIG.loc[:, 'DOWNSTREAM'] = CONFIG['DOWNSTREAM'] + UPPER_LIMIT
CONFIG.loc[:, 'UPSTREAM'] = CONFIG['UPSTREAM'].map('{:,.0f}'.format)
CONFIG.loc[:, 'DOWNSTREAM'] = CONFIG['DOWNSTREAM'].map('{:,.0f}'.format)
CONFIG['NEW_COORD'] = CONFIG['CHR'] + ':' + CONFIG['UPSTREAM'] + '-' + CONFIG['DOWNSTREAM']
CONFIG[['EVENT', 'NEW_COORD', 'COORD']]. \
    to_csv(f'{OUTDIR}/config.tab', sep='\t', index=False, header=None)

# Prepare palette file
pandas.DataFrame([SASHIMI_CONTROL, SASHIMI_CONTRAST]). \
    to_csv(f'{OUTDIR}/palette.txt', index=False, header=None)
