## Quality check of trimmed reads
rule fastqc_after_trimming:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        unpack(get_fastq2)
    output:
        r1 = os.path.join(config["outdir"],"qc", "fastqc", "{sample}.1P_fastqc.zip"),
        r2 = os.path.join(config["outdir"],"qc", "fastqc", "{sample}.2P_fastqc.zip"),
    params:
        outdir = os.path.join(config["outdir"], "qc", "fastqc")
    threads:
        2
    shell:
        """
        fastqc --quiet --outdir {params.outdir} --noextract -f fastq {input} -t {threads}
        """

## Mapping rates etc.
rule bamstats:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        os.path.join(config["outdir"], "bam", "{sample}.bam")
    output:
        os.path.join(config["outdir"], "qc", "bamtools", "{sample}_bamtools.stats")
    threads:
        1
    shell:
        """
        bamtools stats -in {input} | grep -v "*" > {output}
        """

## Make bed windows for window-based metrics: coverage
rule make_ref_window:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        os.path.join(config["outdir"],"ref","ref.fasta.fai")
    output:
        os.path.join(config["outdir"],"ref","ref.window.bed")
    params:
        w_size = config["w_size"]
    threads:
        1
    shell:
        """
        bedtools makewindows -g {input} -w {params.w_size} > {output}
        """

## Coverage
rule genomeCov:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        bam = os.path.join(config["outdir"], "bam", "{sample}.bam"),
        bed = os.path.join(config["outdir"],"ref","ref.window.bed")
    output:
        os.path.join(config["outdir"], "qc", "coverage", "{sample}_coverage.txt")
    threads:
        20
    shell:
        """
        bedtools coverage \
            -a {input.bed} -b {input.bam} \
            > {output}
        """

## vcf stats
rule vcf_stats:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        vcf = os.path.join(config["outdir"],"vcf","all.vcf.gz")
    output:
        vcf_stat = os.path.join(config["outdir"],"qc","bcftools_stats", "vcf.stats")
    threads:
        1
    shell:
        """
        bcftools stats {input.vcf} > {output.vcf_stat}
        """

## Visualize vcf stats
rule plot_vcfstats:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        vcf_stat = os.path.join(config["outdir"],"qc","bcftools_stats","vcf.stats")
    output:
        plot_vcfstats = os.path.join(config["outdir"],"qc","bcftools_stats","plot-vcfstats.log")
    params:
        outdir = os.path.join(config["outdir"],"qc","bcftools_stats")
    threads:
        1
    shell:
        """
        plot-vcfstats \
            -p {params.outdir} \
            --no-PDF \
            {input.vcf_stat}
        """

## Quality check of vcf files
rule qc_vcf:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        vcf = os.path.join(config["outdir"],"vcf","all.vcf.gz")
    output:
        os.path.join(config["outdir"], "qc", "vcftools", "heterozygosity.het")
    params:
        outdir = os.path.join(config["outdir"], "qc", "vcftools")
    threads:
        1
    shell:
        """
        vcftools --gzvcf {input.vcf} --freq2 --out {params.outdir}/allele_frequency --max-alleles 2
        vcftools --gzvcf {input.vcf} --depth --out {params.outdir}/depth_per_indv
        vcftools --gzvcf {input.vcf} --site-mean-depth --out {params.outdir}/depth_per_site
        vcftools --gzvcf {input.vcf} --site-quality --out {params.outdir}/quality_per_site
        vcftools --gzvcf {input.vcf} --missing-indv --out {params.outdir}/missing_rate_per_indv
        vcftools --gzvcf {input.vcf} --missing-site --out {params.outdir}/missing_rate_per_site
        vcftools --gzvcf {input.vcf} --het --out {params.outdir}/heterozygosity
        """

## Summarize all qc files using multiqc
rule multiqc:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        expand(os.path.join(config["outdir"],"qc", "fastqc", "{sample}.{R}_fastqc.zip"), sample= sample_sheet.index, R=["1P", "2P"]),
        expand(os.path.join(config["outdir"], "qc", "bamtools","{sample}_bamtools.stats"), sample= sample_sheet.index),
        expand(os.path.join(config["outdir"],"qc","fastp","{sample}.fastp.json"), sample= sample_sheet.index),
        os.path.join(config["outdir"],"qc","bcftools_stats","vcf.stats")
    output:
        os.path.join(config["outdir"],"qc", "multiqc", "multiqc_report.html")
    params:
        input_dir = os.path.join(config["outdir"], "qc"),
        output_dir = os.path.join(config["outdir"], "qc", "multiqc"),
    threads:
        1
    shell:
        """
        multiqc \
        -o {params.output_dir} \
        {params.input_dir}
        """
