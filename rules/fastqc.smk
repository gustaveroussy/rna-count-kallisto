"""
This rule collects metrics on raw fastq files.
More information at:
https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/fastqc.html
"""
rule fastqc:
    input:
        lambda wildcards: fq_root_dict[wildcards.sample]
    output:
        html = report(
            "qc/fastqc/{sample}_fastqc.html",
            caption="../report/fastqc.rst",
            category="Quality Controls"
        ),
        zip = "qc/fastqc/{sample}_fastqc.zip"
    params:
        ""
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 256, 768)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 20, 120)
        )
    version: "1.0"
    wildcard_constraints:
        sample = r"[^/]+"
    log:
        "logs/fastqc/{sample}.log"
    message:
        "Controling quality of {wildcards.sample} fastq file with FastQC"
    wrapper:
        f"{swv}/bio/fastqc"
