#!/usr/bin/env bash

INPUT="${snakemake_input[0]}"
CONFIG="${snakemake_input[1]}"
PALETTE="${snakemake_input[2]}"
GTF="${snakemake_input[3]}"
OUTDIR="${snakemake_input[4]}"
MIN_COV=${snakemake_params[sashimi_min_cov]}
AGG="${snakemake_params[sashimi_agg_func]}"

while read -r eventID new_coords coords; do
  echo "Plotting $eventID..."
  python /ggsashimi.py \
    --bam "$INPUT" \
    --coordinates "$new_coords" \
    --gtf "$GTF" \
    --out-prefix "${OUTDIR}/sashimi_${eventID}_${coords}" \
    --min-coverage "$MIN_COV" \
    --palette "$PALETTE" \
    --color-factor 3 \
    --overlay 3 \
    --aggr "$AGG" \
    --alpha .6 \
    --fix-y-scale \
    --out-resolution 600 \
    --shrink \
    --ann-height 5 \
    --base-size 12
done <"$CONFIG"
