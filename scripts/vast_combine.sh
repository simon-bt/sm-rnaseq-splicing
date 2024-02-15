#!/usr/bin/env bash

vast-tools combine \
  -sp "${snakemake_params[species]}" \
  --output results/vast_out/ \
  --cores ${snakemake[threads]} \
  "${snakemake_params[opts]}" 2>logs/vt-combine.log
