#!/usr/bin/env bash
set -e

# Edit the sample names, FASTQ paths, references and CPU cores before running.
# FASTQ sample prefixes must match the corresponding filenames.
sample_id="Sample01"
fastq_sample_m="mRNA"
fastq_sample_t="TCR"
fastqs_m="/path/to/mRNA_fastqs"
fastqs_t="/path/to/TCR_fastqs"
reference_m="/path/to/refdata-gex-GRCh38-2020-A"
reference_t="/path/to/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0"
cores=32

# Gene expression. Add --create-bam=true if required by your Cell Ranger version.
cellranger count \
  --id="${sample_id}_mRNA" \
  --fastqs="$fastqs_m" \
  --sample="$fastq_sample_m" \
  --transcriptome="$reference_m" \
  --localcores="$cores"

# T-cell receptor reconstruction.
cellranger vdj \
  --id="${sample_id}_TCR" \
  --fastqs="$fastqs_t" \
  --sample="$fastq_sample_t" \
  --reference="$reference_t" \
  --localcores="$cores" \
  --chain=TR


# Run once after all individual TCR libraries have finished.
# Edit the CSV path; see sampleTcr.example.csv for the required columns.
aggregation_id="TCR"
csv="/path/to/sampleTcr.csv"
cores=32

cellranger aggr \
  --id="$aggregation_id" \
  --csv="$csv" \
  --localcores="$cores"
