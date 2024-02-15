#!/usr/bin/env bash

if [ ${snakemake_params[paired_end]} ]; then
  STAR --readFilesIn "${snakemake_input[0]}" "${snakemake_input[1]}" \
    --genomeDir "${snakemake_input[2]}" \
    --genomeLoad NoSharedMemory \
    --readFilesCommand zcat \
    --twopassMode Basic \
    --outFilterType BySJout \
    --outSAMattributes Standard \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "results/star_out/${snakemake_wildcards[sample]}." \
    --runThreadN ${snakemake[threads]}
else
  STAR --readFilesIn "${snakemake_input[0]}" \
    --genomeDir "${snakemake_input[1]}" \
    --genomeLoad NoSharedMemory \
    --readFilesCommand zcat \
    --twopassMode Basic \
    --outFilterType BySJout \
    --outSAMattributes Standard \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "results/star_out/${snakemake_wildcards[sample]}." \
    --runThreadN ${snakemake[threads]}
fi &&
samtools index --threads ${snakemake[threads]} \
  "results/star_out/${snakemake_wildcards[sample]}.Aligned.sortedByCoord.out.bam" &&
chmod +xr results/star_out/${snakemake_wildcards[sample]}.Aligned.sortedByCoord.out.*
