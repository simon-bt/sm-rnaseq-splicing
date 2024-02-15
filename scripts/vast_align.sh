#!/usr/bin/env bash

if [ ${snakemake_params[paired_end]} ]; then
  vast-tools align "${snakemake_input[0]}" "${snakemake_input[1]}" \
    --sp "${snakemake_params[species]}" \
    --output results/vast_out/ \
    --name ${snakemake_wildcards[sample]} \
    --cores ${snakemake[threads]} 2>"logs/vt-align-${snakemake_wildcards[sample]}.log"
else
  vast-tools align "${snakemake_input[0]}" \
    --sp "${snakemake_params[species]}" \
    --output results/vast_out/ \
    --name ${snakemake_wildcards[sample]} \
    --cores ${snakemake[threads]} 2>"logs/vt-align-${snakemake_wildcards[sample]}.log"
fi
