rule aggregate_kallisto_counts:
    input:
        quants = expand(
            "pseudo_mapping/{sample}",
            sample=fq_root_dict.keys()
        )
