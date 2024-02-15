import pandas
from modules.process import ProcessData

# DO NOT MODIFY - Variables for data visualization
VASTOUT = snakemake.input[0]
DF_ALL = snakemake.input[1]
DF_DIFF = snakemake.input[2]
IMPACT = snakemake.input[3]
MAPPING = snakemake.input[4]

DPSI = snakemake.params['min_dpsi']
SAMPLES_A = snakemake.params['samples_a']
SAMPLES_B = snakemake.params['samples_b']
GROUP_A = snakemake.params['groupA']
GROUP_B = snakemake.params['groupB']
SAMPLES = snakemake.params['samples']
COLOR_GROUPA = snakemake.params['color_groupA']
COLOR_GROUPB = snakemake.params['color_groupB']
COLOR_REG = snakemake.params['color_reg']
COLOR_UNREG = snakemake.params['color_unreg']
OUTDIR = f"results/{GROUP_B}_vs_{GROUP_A}"

# ==========
# Load and prepare data
# ==========

# Load vast-tools table
VASTOUT = pandas.read_csv(VASTOUT, sep="\t")
# Load table for all splicing events passing coverage and default criteria
VAST_ALL = pandas.read_csv(DF_ALL, sep="\t")
# Load table for regulated splicing events passing coverage and dPSI distribution criteria
VAST_DIFF = pandas.read_csv(DF_DIFF, sep="\t")

# Process protein impact dataframe
PROT_IMPACT = pandas. \
    read_csv(IMPACT, sep="\t"). \
    rename(columns={'EventID': 'EVENT'})
prot_mapping = {
    # Protein isoform ontologies
    "Alternative protein isoforms": "ProteinIsoforms",
    "Alternative protein isoforms (Ref)": "ProteinIsoforms",
    "Alternative protein isoforms (No Ref)": "ProteinIsoforms",
    "Alternative protein isoforms (Ref, Alt. Stop)": "ProteinIsoforms",
    "Alternative protein isoforms (No Ref, Alt. Stop)": "ProteinIsoforms",
    "Alternative protein isoforms (Ref, Alt. ATG (>10 exons))": "ProteinIsoforms",
    "Alternative protein isoforms (No Ref, Alt. ATG)": "ProteinIsoforms",
    "Protein isoform when splice site is used (Ref)": "ProteinIsoforms",
    "Protein isoform when splice site is used (No Ref)": "ProteinIsoforms",
    "Protein isoform when splice site is used (No Ref, Alt. Stop)": "ProteinIsoforms",
    "Protein isoform when splice site is used (No Ref, Alt. ATG)": "ProteinIsoforms",
    # ORF-disrupting ontologies
    "ORF disruption upon sequence exclusion": "ORF-disr.",
    "ORF disruption upon sequence inclusion": "ORF-disr.",
    "ORF disruption upon sequence exclusion (Ref, Alt. ATG (<=10 exons))": "ORF-disr.",
    "ORF disruption upon sequence inclusion (Alt. Stop)": "ORF-disr.",
    "ORF disruption when splice site is used (sequence exclusion)": "ORF-disr.",
    "ORF disruption when splice site is used (sequence inclusion)": "ORF-disr.",
    "ORF disruption upon sequence inclusion (1st CDS intron)": "ORF-disr.",
    # Other ontologies
    "5' UTR": "UTRs",
    "3' UTR": "UTRs",
    "NonCoding": "NonCoding",
    "In the CDS, with uncertain impact": "Uncertain"
}
PROT_IMPACT['IMPACT'] = PROT_IMPACT['ONTO'].map(prot_mapping)
PROT_IMPACT.drop('ONTO', axis=1, inplace=True)

# Load and process EventID-GeneID mapping dataframe
EVENT_MAPPING = pandas. \
    read_csv(MAPPING, sep='\t'). \
    rename(columns={'EventID': 'EVENT', 'GeneID': 'ENSEMBL_ID'})

# ==========
# Process vast-tools data
# ==========

# Instantiate the object
ProcessObj = ProcessData(
    vast_out=VASTOUT,
    df_all=VAST_ALL,
    df_diff=VAST_DIFF,
    samples=SAMPLES,
    group_a=GROUP_A,
    group_b=GROUP_B
)

# Get processed dataframes
DATA_ALL = ProcessObj.alt_all.copy()
DATA_ALL_EXIR = ProcessObj.alt_all.copy(). \
    query("TYPE != \'ALT-D\' & TYPE != \'ALT-A\'")
DATA_DIFF = ProcessObj.alt_diff.copy(). \
    merge(PROT_IMPACT, on='EVENT', how='left'). \
    merge(EVENT_MAPPING, on='EVENT', how='left')
DATA_DIFF_EXIR = DATA_DIFF. \
    query("TYPE != \'ALT-D\' & TYPE != \'ALT-A\'")

# ==========
# Visualize alternative splicing data
# ==========

# Plot inclusion data per event type
inclusion_ex = ProcessObj.scatter_quantlevel(
    df_all=DATA_ALL_EXIR,
    event_type="EX",
    group_a=GROUP_A,
    group_b=GROUP_B,
    color_reg=COLOR_REG,
    color_unreg=COLOR_UNREG
)
inclusion_mic = ProcessObj.scatter_quantlevel(
    df_all=DATA_ALL_EXIR,
    event_type="MIC",
    group_a=GROUP_A,
    group_b=GROUP_B,
    color_reg=COLOR_REG,
    color_unreg=COLOR_UNREG
)
inclusion_ri = ProcessObj.scatter_quantlevel(
    df_all=DATA_ALL_EXIR,
    event_type="RI",
    group_a=GROUP_A,
    group_b=GROUP_B,
    color_reg=COLOR_REG,
    color_unreg=COLOR_UNREG
)

# Plot dPSI data
dpsi_dist = ProcessObj.distribution_dpsi(
    df_diff=DATA_DIFF_EXIR,
    dpsi=DPSI,
    group_a=GROUP_A,
    group_b=GROUP_B,
    color_reg=COLOR_REG
)

# Plot inclusion data for regulated events
inclusion_dist = ProcessObj.violin_quantlevel(
    df_all=DATA_ALL_EXIR,
    group_a=GROUP_A,
    group_b=GROUP_B,
    color_group_a=COLOR_GROUPA,
    color_group_b=COLOR_GROUPB
)

# Plot Heat Maps
hmap_standardized = ProcessObj.heatmap(
    df_diff=DATA_DIFF_EXIR,
    samples=SAMPLES,
    zscore=False
)
hmap_zscores = ProcessObj.heatmap(
    df_diff=DATA_DIFF_EXIR,
    samples=SAMPLES,
    zscore=True
)

# Plot ridge dPSI distribution
dpsi_ridge = ProcessObj.ridge_dpsi(
    df_diff=DATA_DIFF_EXIR,
    group_a=GROUP_A,
    group_b=GROUP_B
)

# Plot protein impact data for regulated events
protein_impact = ProcessObj.pct_protimpact(
    df_diff=DATA_DIFF_EXIR
)

# Plot percentage plot of truly regulated events
percentage_dpsi = ProcessObj.pct_dpsi(
    df_all=DATA_ALL_EXIR,
    dpsi=DPSI,
    color_reg=COLOR_REG,
    color_unreg=COLOR_UNREG
)

# Plot quantification levels for regulated events
events_plots = ProcessObj.event_quantlevel(
    df_diff=DATA_DIFF_EXIR,
    samples_a=SAMPLES_A,
    samples_b=SAMPLES_B,
    group_a=GROUP_A,
    group_b=GROUP_B,
    color_group_a=COLOR_GROUPA,
    color_group_b=COLOR_GROUPB
)

# ==========
# Save figures
# ==========

# Plotly figures
figures_plotly = {
    "Scatter_EX": inclusion_ex,
    "Scatter_MIC": inclusion_mic,
    "Scatter_RI": inclusion_ri,
    "DistDPSI": dpsi_dist,
    "LevelsDiff": inclusion_dist,
    "ProteinImpact": protein_impact,
    "PercentageDPSI": percentage_dpsi,
}
for key in figures_plotly:
    figures_plotly[key].write_image(
        file=f"{OUTDIR}/figures/splicing/{key}.pdf",
        format='pdf'
    )
for key in events_plots.keys():
    events_plots[key].write_image(
        file=f"{OUTDIR}/figures/regulated_events/{key}.pdf",
        format='pdf'
    )

# Seaborn figures
figures_seaborn = {
    "RidgeDPSI": dpsi_ridge,
    "HMapLevels_Standardized": hmap_standardized,
    "HMapLevels_ZScores": hmap_zscores,
}
for key in figures_seaborn:
    figures_seaborn[key].savefig(
        f"{OUTDIR}/figures/splicing/{key}.pdf",
        format='pdf',
        dpi=600
    )

# ==========
# Save data
# ==========

# Save processed alternative splicing data
DATA_ALL.to_csv(
    f"{OUTDIR}/data/DataAll-{GROUP_B}-vs-{GROUP_A}_dPSI-{DPSI}.tsv",
    sep="\t",
    index=False)
DATA_DIFF.to_csv(
    f"{OUTDIR}/data/DataDiff-{GROUP_B}-vs-{GROUP_A}_dPSI-{DPSI}.tsv",
    sep="\t",
    index=False)

# Save inclusion data per event type
DFS_ALL = ProcessObj.data_per_event(DATA_ALL)
for key in DFS_ALL.keys():
    df_events = DFS_ALL[key].sort_values("DIFF", ascending=False)
    df_events.to_csv(
        f"{OUTDIR}/data/DataAll-{GROUP_B}-vs-{GROUP_A}_dPSI-{DPSI}_Events-{key}.tsv",
        sep="\t",
        index=False)

# Save dPSI data per event type
DFS_DIFF = ProcessObj.data_per_event(DATA_DIFF)
for key in DFS_DIFF.keys():
    df_events = DFS_DIFF[key]
    df_events.to_csv(
        f"{OUTDIR}/data/DataDiff-{GROUP_B}-vs-{GROUP_A}_dPSI-{DPSI}_Events-{key}.tsv",
        sep="\t",
        index=False)
