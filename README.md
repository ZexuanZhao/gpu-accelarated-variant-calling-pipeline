# Variant Calling from Illumina Reads using GPU

Main project repository is [here](https://github.com/ZexuanZhao/Pegoscapus-hoffmeyeri-sp.A-genome-paper/tree/main).

## Description:
 - A GPU-accelarated snakemake workflow that calls variants from multi-sample illumina reads using Deepvariant and GLnexus

## Files to prepare:
 - A sample sheet - sample_sheet.csv: a comma delimited file with 3 columns (no column name):
   - sample, path to illumina read1, path to illumina read2
 - Modify configuration file - config.yaml:
   - reference:  path to reference fasta file
   - sample_sheet: path to the sample sheet prepared above
   - outdir: path to the output directory
   - suffix: illumina reads' suffix of forward reads and reverse reads. For example:
     - `test1_R1.fastq.gz` and `test1_R2.fastq.gz` should be ["_R1","_R2"]
   - cpu: number of cores provided to the pipeline, should be the same as the command line parameter
   - w_size: non-overlapping window size of reporting average depth along the genome.

## Environment:
 - Make sure snakemake and singularity is installed in current environment.

## Notes:
 - `clara-parabricks` now require 38Gb memory for `fq2bam`. Therefore, `--low-memory` option is used.
## Usage:
`snakemake --use-conda --use-singularity --singularity-args '--nv -B .:/dum' --cores ncpu --resources gpus=ngpu`
