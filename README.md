# Variant Calling from Illumina Reads using GPU

The original project repository is [here](https://github.com/ZexuanZhao/Pegoscapus-hoffmeyeri-sp.A-genome-paper/tree/main).

## Description:
 - A GPU-accelarated snakemake workflow that calls variants from multi-sample illumina reads using Deepvariant and GLnexus

## Files to prepare:
 - A sample sheet - sample_sheet.csv: a comma delimited file with 3 columns (no column name):
   - sample, path to illumina read1, path to illumina read2
 - Modify configuration file - `configuration/config.yaml`:
   - project: a name for project
   - reference:  path to reference fasta file
   - sample_sheet: path to the sample sheet prepared above
   - outdir: path to the output directory
   - clara-parabricks: image path to clara-parabricks, e.g. "docker://nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1"
   - deepvariant: image path to deepvariant, e.g. "docker://google/deepvariant:1.8.0_saved_model-gpu"
   - w_size: non-overlapping window size of reporting average depth along the genome.

## Environment:
 - Make sure snakemake and singularity is installed in current environment.

## Notes:
 - `clara-parabricks` now require 38Gb memory for `fq2bam`. Therefore, `--low-memory` option is used.
 - In snakemake commandline, [ncpu] should be larger than 20 as all resource usages are hardcoded.
## Usage:
`snakemake --use-conda --use-singularity --singularity-args '--nv -B .:/dum' --cores [ncpu] --resources gpus=[ngpu]`
