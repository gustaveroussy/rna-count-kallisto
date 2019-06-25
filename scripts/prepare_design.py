#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""
This script aims to prepare the list of files to be processed
by the rna-count-kallisto pipeline
"""

import pandas as pd                   # Parse TSV files
from argparse import ArgumentParser   # Parse command line
from pathlib import Path              # Paths related methods
import sys                            # System related methods
from typing import Generator          # Type hints


def search_fq(fq_dir: Path,
              recursive: bool = False) -> Generator[str, str, None]:
    """
    Iterate over a directory and search for fastq files
    """
    for path in fq_dir.iterdir():
        if path.is_dir():
            if recursive is True:
                yield from search_fq(path, recursive)
            else:
                continue

        if path.name.endswith(("fq", "fq.gz", "fastq", "fastq.gz")):
            yield path


if __name__ == '__main__':
    main_parser = ArgumentParser(
        description="Prepare the design file for rna-count-kallisto",
        epilog="Each tool belong to their respective authors"
    )

    main_parser.add_argument(
        "path",
        help="Path to the directory containing fastq files",
        type=str
    )

    main_parser.add_argument(
        "-s", "--single",
        help="The samples are single ended rnaseq reads, not pair ended",
        action="store_true"
    )

    main_parser.add_argument(
        "-r", "--recursive",
        help="Recursively search in sub-directories for fastq files",
        action="store_true"
    )

    main_parser.add_argument(
        "-o", "--output",
        help="Path to output file (default: %(default)s)",
        type=str,
        default="design.tsv"
    )

    args = main_parser.parse_args()

    fq_files = sorted(list(search_fq(Path(args.path), args.recursive)))
    fq_dict = {}
    if args.single is True:
        for fq in fq_files:
            fq_dict[fq.name] = {
                "Sample_id": fq.stem,
                "Upstream_file": fq.absolute()
            }
    else:
        for fq1, fq2 in zip(fq_files[0::2], fq_files[1::2]):
            fq_dict[fq1.name] = {
                "Sample_id": fq1.stem,
                "Upstream_file": fq1.absolute(),
                "Downstream_file": fq2.absolute()
            }

    data = pd.DataFrame(fq_dict).T
    print(data.head(), file=sys.stderr)
    data.to_csv(args.output, sep="\t", index=False)
