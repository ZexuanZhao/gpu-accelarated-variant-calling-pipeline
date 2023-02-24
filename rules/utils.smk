def get_fastq(wildcards):
    fastqs = sample_sheet.loc[(wildcards), ["r1", "r2"]].dropna()
    return {"r1": fastqs.r1, "r2": fastqs.r2}

def get_fastq2(wildcards):
    return {"r1": os.path.join(config["outdir"], "trimmed_reads", wildcards.sample + ".1P.fastq.gz"),
            "r2": os.path.join(config["outdir"], "trimmed_reads", wildcards.sample + ".2P.fastq.gz")}

def get_fastq2_basename(wildcards):
    return {"r1_name": wildcards.sample + ".1P.fastq.gz",
            "r2_name": wildcards.sample + ".2P.fastq.gz"}