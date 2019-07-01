"""
This rule runs MultiQC in order to collect metrics on most of our tools and
raw files: Fastq + Kallisto. We need to include the fasta reference for
the report option only.
More information at:
https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/multiqc.html
"""
rule multiqc:
    input:
        # Include FastQC reports
        html = expand(
            "qc/fastqc/{sample}_fastqc.{ext}",
            sample=fq_root_dict.keys(),
            ext=["html", "zip"]
        ),
        # Include Kallisto quants
        quants = expand(
            "pseudo_mapping/{sample}",
            sample=sample_id_list
        )
    output:
        report(
            "qc/multiqc_report.html",
            caption="../report/multiqc.rst",
            category="Quality Controls"
        )
    params: ""
    threads: 1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 256, 768)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 10, 60)
        )
    version: "1.0"
    log:
        "logs/multiqc.log"
    message:
        "Gathering quality reports with MultiQC"
    wrapper:
        f"{swv}/bio/multiqc"
