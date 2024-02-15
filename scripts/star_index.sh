#!/usr/bin/env bash

STAR --runMode genomeGenerate \
  --genomeDir "results/star_out/star_index" \
  --genomeFastaFiles "${snakemake_input[0]}" \
  --sjdbGTFfile "${snakemake_input[1]}" \
  --sjdbOverhang "${snakemake_params[read_length]}" \
  --runThreadN ${snakemake[threads]}
