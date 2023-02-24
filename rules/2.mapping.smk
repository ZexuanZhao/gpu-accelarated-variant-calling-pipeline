rule copy_reference:
    input:
        original_ref = config["reference"]
    output:
        copied_ref = os.path.join(config["outdir"], "ref", "ref.fasta")
    shell:
        """
        cp {input.original_ref} {output.copied_ref}
        """

rule index_reference:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        ref = os.path.join(config["outdir"], "ref", "ref.fasta")
    output:
        ref_bwt = os.path.join(config["outdir"], "ref", "ref.fasta.bwt")
    shell:
        """
        bwa index {input.ref}
        """
rule bwa_gpu:
    input:
        unpack(get_fastq2),
        ref_bwt = os.path.join(config["outdir"], "ref", "ref.fasta.bwt"),
    output:
        os.path.join(config["outdir"], "bam", "{sample}.bam")
    params:
        reads_dict = get_fastq2_basename,
        ref_path = os.path.abspath(os.path.join(config["outdir"], "ref")),
        read_path = os.path.abspath(os.path.join(config["outdir"], "trimmed_reads")),
        tmp_path = os.path.abspath(os.path.join(config["outdir"], "tmp", "bam", "{sample}")),
        out_bam_path = os.path.join(config["outdir"], "bam"),
        out_bam_name = "{sample}.bam"
    log:
        os.path.join(config["outdir"],"logs","mapping","{sample}.log")
    threads:
        config["cpu"]
    shell:
        """
        docker run \
            --gpus all \
            -w /workdir \
            --volume {params.ref_path}:/ref_dir \
            --volume {params.read_path}:/read_dir\
            --volume {params.tmp_path}:/outputdir \
            nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1 \
            pbrun fq2bam \
                --ref /ref_dir/ref.fasta \
                --in-fq /read_dir/{params.reads_dict[r1_name]} /read_dir/{params.reads_dict[r2_name]} \
                --out-bam /outputdir/{params.out_bam_name} \
            > {log} \
            2>{log}
        cp {params.tmp_path}/* {params.out_bam_path}
        """
