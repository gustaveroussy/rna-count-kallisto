"""
While other .smk files contains rules and pure snakemake instructions, this
one gathers all the python instructions surch as config mappings or input
validations.
"""

from snakemake.utils import validate
from typing import Any, Dict, List

import os.path as op    # Path and file system manipulation
import os               # OS related operations
import pandas as pd     # Deal with TSV files (design)
import sys              # System related operations

# Snakemake-Wrappers version
swv = "0.35.1"
# github prefix
git = "https://bitbucket.org/tdayris/snakemake-wrappers/raw"

# Loading configuration
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

# Loading deisgn file
design = pd.read_csv(
    config["design"],
    sep="\t",
    header=0,
    index_col=None,
    dtype=str
)
design.set_index(design["Sample_id"])
validate(design, schema="../schemas/design.schema.yaml")

report: "../report/general.rst"


def fq_link() -> Dict[str, str]:
    """
    This function takes the "samples" described in config and returns
    a dictionnary with:
    sample file name : sample path
    """
    # Will cause KeyError on single stranded RNA-Seq analysis
    # Better ask forgiveness than permission !
    try:
        # Paired-ended case
        fq_list = chain(design["Upstream_file"], design["Downstream_file"])
    except KeyError:
        # Single ended case
        fq_list = design["Upstream_file"]
    finally:
        return {
            op.basename(fq): op.realpath(fq)
            for fq in fq_list
        }


def fq_root() -> Dict[str, str]:
    """
    This function takes the fastq file list and returns the root
    name corresponding to a fastq file
    sample name: sample link path
    """
    # For now, bz2 compression is not taken into account.
    possible_ext = ("fq", "fastq", "fq.gz", "fastq.gz")

    # Will cause KeyError on single stranded RNA-Seq analysis
    # Better ask forgiveness than permission !
    try:
        # Paired-ended case
        fq_list = chain(design["Upstream_file"], design["Downstream_file"])
    except KeyError:
        # Single ended case
        fq_list = design["Upstream_file"]

    # Build final result
    result = {}
    for fq in fq_list:
        # I always love writing these crazy for-break-else!
        for ext in possible_ext:
            if fq.endswith(ext):
                # Extension removal
                base = op.basename(fq)[:-(len(ext) + 1)]
                result[base] = f"raw_data/{op.basename(fq)}"
                break
        else:
            raise ValueError(f"Could not remove ext: {fq}")

    return result


def ref_link() -> Dict[str, str]:
    """
    This function takes the "ref" described in config and returns
    a dictionnary with:
    ref file name : ref path
    """
    # If not GTF is provided, error will be raised.
    try:
        # Case GTF is provided
        fasta = config["ref"]["fasta"]
        gtf = config["ref"]["gtf"]
        return {
            op.basename(fasta): op.realpath(fasta),
            op.basename(gtf): op.realpath(gtf)
        }
    except KeyError:
        # Case GTF is missing
        fasta = config["ref"]["fasta"]
        return {
            op.basename(fasta): op.realpath(fasta)
        }
    except TypeError:
        # Case GTF is missing
        fasta = config["ref"]["fasta"]
        return {
            op.basename(fasta): op.realpath(fasta)
        }


def fq_pairs() -> Dict[str, str]:
    """
    This function returns a sample ID and
    the corresponding fastq files.
    """
    # Will cause KeyError on single stranded RNA-Seq analysis
    # Better ask forgiveness than permission !
    try:
        # Paired end case
        iterator = zip(
            design["Sample_id"],
            design["Upstream_file"],
            design["Downstream_file"]
        )
        return {
            name: [
                f"raw_data/{op.basename(fq1)}",
                f"raw_data/{op.basename(fq2)}"
            ]
            for name, fq1, fq2 in iterator
        }
    except KeyError:
        # Single end case
        iterator = zip(
            design["Sample_id"],
            design["Upstream_file"]
        )
        raise
        return {
            name: [f"raw_data/{op.basename(fq1)}"]
            for name, fq1 in iterator
        }


def refs_pack() -> Dict[str, str]:
    """
    Return a dictionnary with references
    """
    # Will cause KeyError if no GTF is given.
    # Better ask forgiveness than permission !
    try:
        # GTF is present
        return {
            "fasta": f"genome/{op.basename(config['ref']['fasta'])}",
            "gtf": f"genome/{op.basename(config['ref']['gtf'])}"
        }
    except KeyError:
        # No GTF provided !
        return {
            "fasta": f"genome/{op.basename(config['ref']['fasta'])}"
        }
    except TypeError:
        # No GTF provided !
        return {
            "fasta": f"genome/{op.basename(config['ref']['fasta'])}"
        }


def fq_pairs_w(wildcards) -> Dict[str, str]:
    """
    Dynamic wildcards call for snakemake.
    """
    return {"fastq": fq_pairs_dict[wildcards.sample]}


def sample_id() -> List[str]:
    """
    Return the list of samples identifiers
    """
    return design["Sample_id"].tolist()


def kallisto_quant_extra() -> str:
    """
    Return the corrected list of parameters for kallist quant
    """
    try:
        return (f"{config['params']['kallisto_quant_extra']} "
                f"--gtf {str(config['gtf'])}")
    except KeyError:
        return f"{config['params']['kallisto_quant_extra']}"


def get_targets() -> Dict[str, Any]:
    """
    This function returns the targets of Snakemake
    following the requests from the user.
    """
    targets = {}
    if config["workflow"]["fastqc"] is True:
        targets["fastqc"] = expand(
            "qc/fastqc/{samples}_fastqc.{ext}",
            samples=fq_root_dict.keys(),
            ext=["html", "zip"]
        )
    if config["workflow"]["multiqc"] is True:
        targets["multiqc"] = "qc/multiqc_report.html"
    if config["workflow"]["aggregate"] is True:
        targets["aggregation"] = [
            "aggregated_kallisto_counts/merged_est_counts.tsv",
            "aggregated_kallisto_counts/merged_tpm.tsv"
        ]
    targets["quant"] = expand(
        "pseudo_mapping/{sample}",
        sample=sample_id_list
    )
    print(targets, file=sys.stderr)
    return targets


# We will use these functions multiple times. On large input datasets,
# pre-computing all of these makes Snakemake faster.
fq_link_dict = fq_link()
fq_root_dict = fq_root()
ref_link_dict = ref_link()
fq_pairs_dict = fq_pairs()
refs_pack_dict = refs_pack()
sample_id_list = sample_id()
targets_dict = get_targets()
