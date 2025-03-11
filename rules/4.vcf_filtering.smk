rule filter_biallelic_snp:
    conda:
        os.path.join(workflow.basedir,"envs/envs.yaml")
    input:
        os.path.join(config["outdir"], "vcf", config["project"]+".vcf.gz")
    output:
        os.path.join(config["outdir"], "vcf", config["project"]+".biallelic.vcf.gz")
    shell:
        """
        bcftools view \
            -m2 -M2 \
            -v snps \
            -O z \
            -o {output} \
            {input}
        """
