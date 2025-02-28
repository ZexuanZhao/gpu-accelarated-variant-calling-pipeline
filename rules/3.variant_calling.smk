rule deepvariant:
    input:
        ref = os.path.join(config["outdir"], "ref", "ref.fasta"),
        bam = os.path.join(config["outdir"],"bam","{sample}.bam")
    output:
        vcf = os.path.join(config["outdir"],"vcf", "{sample}.vcf"),
        gvcf = os.path.join(config["outdir"],"vcf","{sample}.gvcf"),
        report = os.path.join(config["outdir"], "qc", "deepvariant", "{sample}.visual_report.html")
    params:
        sample_name = "{sample}",
        ref_path = os.path.abspath(os.path.join(config["outdir"],"ref")),
        bam_path = os.path.abspath(os.path.join(config["outdir"], "bam")),
        bam_name = "{sample}.bam",
        tmp_path = os.path.abspath(os.path.join(config["outdir"], "tmp", "deepvariant", "{sample}")),
        out_vcf_path = os.path.join(config["outdir"],"vcf"),
        out_vcf_name = "{sample}.vcf",
        out_gvcf_name = "{sample}.gvcf",
        out_report_name = "{sample}.visual_report.html"
    log:
        os.path.join(config["outdir"],"logs","deepvariant","{sample}.log")
    threads:
        config["cpu"]
    shell:
        '''
        docker run \
            --rm \
            --gpus all \
            -w /workdir \
            --volume {params.ref_path}:/ref_dir \
            --volume {params.bam_path}:/bam_dir \
            --volume {params.tmp_path}:/outputdir \
            google/deepvariant:1.8.0_saved_model-gpu \
            /opt/deepvariant/bin/run_deepvariant \
                --model_type=WGS \
                --ref=/ref_dir/ref.fasta \
                --reads=/bam_dir/{params.bam_name} \
                --sample_name={params.sample_name} \
                --output_vcf=/outputdir/{params.out_vcf_name} \
                --output_gvcf=/outputdir/{params.out_gvcf_name} \
                --num_shards={threads} \
            > {log} \
            2> {log}
        cp {params.tmp_path}/* {params.out_vcf_path}
        mv {params.out_vcf_path}/{params.out_report_name} {output.report}
        '''

rule gvcf2vcf_deepvariant:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        gvcf = expand(os.path.join(config["outdir"],"vcf","{sample}.gvcf"), sample = sample_sheet.index)
    output:
        vcf = os.path.join(config["outdir"],"vcf", "all.vcf.gz")
    log:
        os.path.join(config["outdir"],"logs","deepvariant","glnexus.log")
    threads:
        config["cpu"]
    shell:
        '''
        rm -rf ./GLnexus.DB
        [ -f $CONDA_PREFIX/lib/libjemalloc.so ] && export LD_PRELOAD=$CONDA_PREFIX/lib/libjemalloc.so
        glnexus_cli \
            --config DeepVariantWGS \
            --threads {threads} \
            {input.gvcf} | \
        bcftools view - | \
        bgzip -c \
        > {output.vcf} \
        2> {log}
        rm -rf ./GLnexus.DB
        '''
