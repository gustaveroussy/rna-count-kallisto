#!/usr/bin/python3.7
# conding: utf-8

import os.path as op
import pandas as pd
import sys

from snakemake.utils import makedirs
from typing import List


def read_kallisto(path: str) -> pd.DataFrame:
    """
    Reader for kallisto files
    """
    return pd.read_csv(
        op.join(path, "abundance.tsv"),
        sep="\t",
        index_col=0,
        header=0,
        dtype={
            0: str,
            1: pd.np.float,
            2: pd.np.float,
            3: pd.np.float,
            4: pd.np.float
        }
    )


def extract_field(*paths: List[str],
                  prefix: str = "",
                  column: str = "TPM") -> pd.DataFrame:
    """
    Return a dataframe containing the required field
    """
    merged_frame = None
    for path in paths:
        print(f"Working on {path}", file=sys.stderr)
        data = read_kallisto(path)
        sample_id = str(path)
        sample_id = sample_id if prefix == "" else sample_id[len(prefix):]

        data = data[[column]]
        data.columns = [sample_id]

        try:
            merged_frame = pd.merge(
                merged_frame,
                data,
                left_index=True,
                right_index=True
            )
        except TypeError:
            merged_frame = data

    merged_frame.fillna(0)
    return merged_frame


if __name__ == '__main__':
    makedirs(snakemake.output)

    for column in ["est_counts", "tpm"]:
        data = extract_field(
            *snakemake.input["quants"],
            prefix="pseudo_mapping/",
            column=column
        )

        print(data.head(), file=sys.stderr)

        data.to_csv(snakemake.output[column], sep="\t")
