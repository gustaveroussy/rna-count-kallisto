"""
This rule computes the kallisto index from a fasta-formatted
transcriptome sequence (no genome)

more information at:
https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/kallisto/index.html
"""
rule kallisto_index:
    input:
        **refs_pack_dict
    output:
        index = "pseudo_mapping/genome_index"
    message:
        "Indexing {input.fasta} with Kallisto"
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 2048 + 10240, 35580)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 15 + 15, 180)
        )
    version: "1.0"
    threads:
        min(config["threads"], 12)
    params:
        extra = config["params"]["kallisto_index_extra"]
    log:
        "logs/kallisto/index.log"
    wrapper:
        f"{swv}/bio/kallisto/index"

"""
This rule performs the quantification step with pseudo-mapping
fromfastq-formatted rnaseq reads and kallisto genome index.

More information at:
https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/kallisto/quant.html
"""
rule kallisto_quant:
    input:
        unpack(fq_pairs_w),
        index = "pseudo_mapping/genome_index"
    output:
        report(
            directory("pseudo_mapping/{sample}"),
            caption="../report/counts.rst",
            category="Counts"
        )
    message:
        "Quantifying {wildcards.sample} with Kallisto"
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 5120 + 2048, 20480)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 15 + 15, 180)
        )
    version: "1.0"
    threads:
        min(config["threads"], 12)
    params:
        extra = kallisto_quant_extra()
    log:
        "logs/kallisto/quant_{sample}.log"
    wrapper:
        f"{swv}/bio/kallisto/quant"
