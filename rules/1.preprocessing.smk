rule trim_reads:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        unpack(get_fastq)
    output:
        r1 = os.path.join(config["outdir"], "trimmed_reads", "{sample}.1P.fastq.gz"),
        r2 = os.path.join(config["outdir"], "trimmed_reads", "{sample}.2P.fastq.gz"),
        r1_unpaired= os.path.join(config["outdir"], "trimmed_reads", "{sample}.1U.fastq.gz"),
        r2_unpaired= os.path.join(config["outdir"], "trimmed_reads", "{sample}.2U.fastq.gz"),
        json_report = os.path.join(config["outdir"], "qc", "fastp", "{sample}.fastp.json"),
        html_report = os.path.join(config["outdir"], "qc", "fastp", "{sample}.fastp.html")
    log:
        os.path.join(config["outdir"],"logs", "fastp", "{sample}.log")
    threads:
        2
    shell:
        """
        fastp \
            --in1 {input.r1} \
            --out1 {output.r1} \
            --in2 {input.r2} \
            --out2 {output.r2} \
            --unpaired1 {output.r1_unpaired} \
            --unpaired2 {output.r2_unpaired} \
            --thread {threads} \
            --json {output.json_report} \
            --html {output.html_report} \
            > {log} \
            2> {log}
        """