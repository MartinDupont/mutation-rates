import argparse
import numpy as np
from typing import Optional

import pandas as pd

from common import parse_fasta_file, loop_over_ncbi_folders


def read_fasta_files(base_directory, lim: Optional[int]=100):
    all_sequences = []

    for file_path, folder_name in loop_over_ncbi_folders(base_directory, lim):
            print(f'reading CDS from {folder_name}')
            with open(file_path, 'r') as file:
                parsed = parse_fasta_file(file, folder_name)
                all_sequences += parsed

    df_sequences = pd.DataFrame(all_sequences)
    return df_sequences


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate daily zspreads")
    parser.add_argument(
        "--dir",
        "-d",
        dest="base_directory",
        help="Where the fasta files are",
        type=str,
    )
    parser.add_argument(
        "--limit",
        "-l",
        dest="limit",
        help="The file where the reference genome is stored",
        type=int,
        default=np.inf
    )

    args = parser.parse_args()
    print(args)

    df = read_fasta_files(args.base_directory, lim=args.limit)
    print(df)