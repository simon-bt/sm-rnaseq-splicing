#!/usr/bin/env bash

curl "${snakemake_config[genome]}" -o "${snakemake_output[0]}.gz"
curl "${snakemake_config[annotation]}" -o "${snakemake_output[1]}.gz"
curl "${snakemake_config[prot_impact]}" -o "${snakemake_output[2]}.gz"
curl "${snakemake_config[event_mapping]}" -o "${snakemake_output[3]}.gz"
curl "${snakemake_config[event_info]}" -o "${snakemake_output[4]}.gz"
gunzip "${snakemake_output[0]}.gz" "${snakemake_output[1]}.gz" "${snakemake_output[2]}.gz" \
  "${snakemake_output[3]}.gz" "${snakemake_output[4]}.gz"
