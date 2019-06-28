rule aggregate_kallisto_counts:
    input:
        quants = expand(
            "pseudo_mapping/{sample}",
            sample=sample_id_list
        )
    output:
        report(
            directory("aggregated_kallisto_counts"),
            caption="../report/aggregation.rst",
            category="Aggregation"
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
