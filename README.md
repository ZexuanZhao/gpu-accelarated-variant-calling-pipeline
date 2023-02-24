Description:
 - A gpu-accelerated variant calling pipeline with quality controls using `Deepvariant` and `GLnexus`.


Files to prepare:
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

Environment:
 - Make sure snakemake is installed in current environment.
 - Docker is required.
 - Install docker image: `nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1`
 - Install docker image: `google/deepvariant:1.4.0-gpu`

To run: `snakemake --cores [cpu] --use-conda`