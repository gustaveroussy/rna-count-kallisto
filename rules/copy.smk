"""
On most clusters, cold and hot storage coexist. Non-expert users might
try to run IO intensive processes on data through cold storage and break
either the pipeline or the mounting points on a cluster. This rule
copies the fastq files.
More information at:
https://github.com/tdayris/yawn/tree/master/SnakemakeWrappers/cp/8.25
"""
rule copy_fastq:
    input:
        lambda wildcards: fq_link_dict[wildcards.files]
    output:
        temp("raw_data/{files}")
    message:
        "Copying {wildcards.files} for further process"
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 128, 512)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 1440, 2832)
        )
    version: "1.0"
    log:
        "logs/copy_{files}.log"
    wildcard_constraints:
        files = r"[^/]+"
    threads: 1
    priority: 1
    params:
        extra = config["params"].get("copy_extra", ""),
        cold_storage = config.get("cold_storage", "NONE")
    wrapper:
        f"{git}/cp/bio/cp"

"""
Same remarks as the above. Here, we copy the reference files.
"""
rule copy_extra:
    input:
        lambda wildcards: ref_link_dict[wildcards.files]
    output:
        temp("genome/{files}")
    message:
        "Copying {wildcards.files} as reference"
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 128, 512)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 1440, 2832)
        )
    version: "1.0"
    log:
        "logs/copy_{files}.log"
    wildcard_constraints:
        files = r"[^/]+"
    threads: 1
    priority: 1
    params:
        extra = config["params"].get("copy_extra", ""),
        cold_storage = config.get("cold_storage", "NONE")
    wrapper:
        f"{git}/cp/bio/cp"
