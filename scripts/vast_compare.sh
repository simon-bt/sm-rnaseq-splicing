#!/usr/bin/env bash

if [ ${snakemake_params[paired]} ]; then
    vast-tools compare "${snakemake_input[0]}" \
    --min_dpsi ${snakemake_params[min_dpsi]} \
    --min_range ${snakemake_params[min_range]} \
    --samplesA "${snakemake_params[samplesA]}" \
    --samplesB "${snakemake_params[samplesB]}" \
    -name_A "${snakemake_params[groupA]}" \
    -name_B "${snakemake_params[groupB]}" \
    --print_dPSI \
    --print_all_ev \
    --outFile "${snakemake_params[groupB]}-vs-${snakemake_params[groupA]}_dPSI-${snakemake_params[min_dpsi]}.tab"
else
    vast-tools compare "${snakemake_input[0]}" \
    --min_dpsi ${snakemake_params[min_dpsi]} \
    --min_range ${snakemake_params[min_range]} \
    --samplesA "${snakemake_params[samplesA]}" \
    --samplesB "${snakemake_params[samplesB]}" \
    -name_A "${snakemake_params[groupA]}" \
    -name_B "${snakemake_params[groupB]}" \
    --print_dPSI \
    --print_all_ev \
    --paired \
    --outFile "${snakemake_params[groupB]}-vs-${snakemake_params[groupA]}_dPSI-${snakemake_params[min_dpsi]}.tab"
fi
