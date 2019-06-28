rule aggregate_kallisto_counts:
    input:
        quants = expand(
            "pseudo_mapping/{sample}",
            sample=sample_id_list
        )
    output:
        est_counts = report(
            "aggregated_kallisto_counts/merged_est_counts.tsv",
            caption="../report/raw_counts.rst",
            category="Aggregated Counts"
        ),
        tpm = report(
            "aggregated_kallisto_counts/merged_tpm.tsv",
            caption="../report/tpm.rst",
            category="Aggregated Counts"
        )
    message:
        "Aggregating all kallisto abundancies"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(
                attempt * len(sample_id_list) * 250, 10240
            )
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 10, 60)
        )
    log:
        "logs/aggregation.log"
    conda:
        "../envs/py37.yaml"
    script:
        "../scripts/aggregate_samples.py"
