import os
import pandas as pd
import pdb
import re
from typing import List, Optional


def parse_fasta_file(file_handle, sample_name) -> List[dict]:
    records = []

    lines = file_handle.readlines()
    for line in lines:
        if not line.startswith(">"):
            continue
        description = line

        #regex = r'\[(.*?)\]'
        regex = r'\[([^][]*(?:\[[^][]*\])*[^][]*)\]'
        matches = re.findall(regex, description)
        record_data = {}

        for match in matches:
            if '=' in match:
                key, value = match.split('=', 1)
                key = key.strip()
                value = value.strip()
                record_data[key] = value


        record_data['sample'] = sample_name
        records.append(record_data)

    return records



def read_fasta_files(base_directory, lim=100):
    all_sequences = []

    count = 0
    for folder_name in os.listdir(base_directory):
        if count > lim:
            break

        folder_path = os.path.join(base_directory, folder_name)
        if os.path.isdir(folder_path):
            file_path = os.path.join(folder_path, "cds_from_genomic.fna")

            if os.path.isfile(file_path):
                print(f'reading CDS from {folder_path}')
                with open(file_path, 'r') as file:
                    parsed = parse_fasta_file(file, folder_name)
                    all_sequences += parsed
                count += 1

    df_sequences = pd.DataFrame(all_sequences)
    return df_sequences


if __name__ == '__main__':
    base_directory = '/Users/martin/Documents/data/ncbi_new/ncbi_dataset/ncbi_dataset/data'
    df = read_fasta_files(base_directory)
    print(df)