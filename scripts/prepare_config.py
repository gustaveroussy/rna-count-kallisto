#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""
This script aims to prepare the configuration file used
by the rna-count-kallisto pipeline
"""

import yaml                           # Parse Yaml files
from argparse import ArgumentParser   # Parse command line
from pathlib import Path              # Paths related methods
import sys                            # System related methods


if __name__ == '__main__':
    main_parser = ArgumentParser(
        description="Prepare config.yaml file for your Snakefile",
        epilog="This script does not make any magic. Please check the prepared"
               " configuration file!",
    )

    main_parser.add_argument(
        "fasta",
        help="Path to the fasta-formatted transcriptome sequence",
        type=str,
    )

    main_parser.add_argument(
        "--design",
        help="Path to design file",
        type=str,
        metavar="PATH",
        default="design.tsv"
    )

    main_parser.add_argument(
        "--workdir",
        help="Path to working directory",
        type=str,
        metavar="PATH",
        default="."
    )

    main_parser.add_argument(
        "--threads",
        help="Maximum number of threads used",
        type=int,
        default=1
    )

    main_parser.add_argument(
        "--singularity",
        help="Docker/Singularity image",
        type=str,
        default="docker://continuumio/miniconda3:4.4.10"
    )

    main_parser.add_argument(
        "--cold_storage",
        help="Space separated list of absolute path to "
             "cold storage mount points",
        nargs="+",
        type=str,
        default=[" "]
    )

    main_parser.add_argument(
        "--gtf",
        help="Path to GTF-formatted genome annotation",
        default=None,
        type=str
    )

    main_parser.add_argument(
        "--no_quality_control",
        help="Do not perform any additional quality controls",
        action="store_true"
    )

    main_parser.add_argument(
        "--aggregate",
        help="Perform sample count aggregation",
        action="store_true"
    )

    main_parser.add_argument(
        "--kallisto_index_extra",
        help="Extra parameters for kallisto index step",
        type=str,
        default="--make-unique"
    )

    main_parser.add_argument(
        "--kallisto_quant_extra",
        help="Extra parameters for kallisto quantification step",
        type=str,
        default="--bias --bootstrap-samples=100"
    )

    args = main_parser.parse_args()
    config_params = {
        "design": args.design,
        "workdir": args.workdir,
        "threads": args.threads,
        "singularity_docker_image": args.singularity,
        "cold_storage": args.cold_storage,
        "ref": {
            "fasta": args.fasta,
            "gtf": args.gtf
        },
        "workflow": {
            "fastqc": not args.no_quality_control,
            "multiqc": not args.no_quality_control,
            "aggregate": args.aggregate
        },
        "params": {
            "kallisto_index_extra": args.kallisto_index_extra,
            "kallisto_quant_extra": args.kallisto_quant_extra
        }
    }

    output_path = Path(args.workdir) / "config.yaml"
    with output_path.open("w") as config_yaml:
        config_yaml.write(yaml.dump(config_params, default_flow_style=False))
