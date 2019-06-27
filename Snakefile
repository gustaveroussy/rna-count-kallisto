import snakemake.utils  # Load snakemake API
import sys              # System related operations

# Python 3.7 is required
if sys.version_info < (3, 7):
    raise SystemError("Please use Python 3.7 or later.")

# Snakemake 5.4.2 at least is required
snakemake.utils.min_version("5.4.2")

include: "rules/common.smk"
include: "rules/copy.smk"
include: "rules/fastqc.smk"
include: "rules/multiqc.smk"
include: "rules/kallisto.smk"
include: "rules/aggregation.smk"

workdir: config["workdir"]
singularity: config["singularity_docker_image"]
localrules: copy_fastq, copy_extra

rule all:
    input:
        **targets_dict
    message:
        "Finishing the Kallisto RNA-Seq quantification pipeline"
