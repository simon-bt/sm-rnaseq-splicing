#!/usr/bin/env bash

singularity build "${snakemake_output[0]}" docker://guigolab/ggsashimi:latest