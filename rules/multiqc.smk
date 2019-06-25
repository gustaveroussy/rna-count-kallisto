"""
This rule runs MultiQC in order to collect metrics on most of our tools and
raw files: Fastq + STAR + Samtools. We need to include the fasta reference for
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
        ),
        # Include fasta ref for reporting.
        fasta = f"genome/{op.basename(config['ref']['fasta'])}"
    output:
        report("qc/report.html",
               caption="../report/quality.rst",
               category="Quality Controls & Mapping")
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
