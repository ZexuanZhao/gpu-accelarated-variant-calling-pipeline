#!/usr/bin/env python
import os.path
import pandas as pd

configfile: "configuration/config.yaml"
sample_sheet = pd.read_csv(config["sample_sheet"],
    dtype=str,
    names = ["sample", "r1", "r2"]).set_index("sample")

wildcard_constraints:
    sample = "|".join(sample_sheet.index)

include: "rules/utils.smk"
include: "rules/0.qc.smk"
include: "rules/1.preprocessing.smk"
include: "rules/2.mapping.smk"
include: "rules/3.variant_calling.smk"
include: "rules/4.vcf_filtering.smk"

rule all:
    input:
        multiqc_report = os.path.join(config["outdir"],"qc","multiqc", config["project"]+"_multiqc_report.html"),
        bam_report = expand(os.path.join(config["outdir"],"qc","coverage","{sample}_coverage.txt"), sample = sample_sheet.index),
        vcf_plot = os.path.join(config["outdir"],"qc","bcftools_stats","plot-vcfstats.log"),
        biallelic_snp = os.path.join(config["outdir"], "vcf", config["project"]+".biallelic.vcf.gz"),
        het = os.path.join(config["outdir"], "qc", "vcftools", config["project"]+".het")
    shell:
        """
        echo "Job done!"
        echo "Use the following command to clean up temporary files (needs sudo):"
        """
