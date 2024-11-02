import os
import pandas as pd
import pdb
import re
from typing import List, Optional

COLUMNS_OF_INTEREST = ['gene', 'locus_tag', 'protein', 'protein_id']

ERROR_COLUMNS = [
    'pseudo',
    'partial',
    'exception',
    'transl_except',
]

def iterate_fasta(file_handle):
    lines = file_handle.readlines()
    if not lines[0].startswith(">"):
        raise RuntimeError("Fasta file not formatted correctly. Doesn't start with >")

    sequence = None
    descriptions = []
    sequences = []
    for line in lines:
        if not line.startswith(">"):
            sequence += line.strip()
        else:
            if sequence is not None:
                sequences.append(sequence)
            sequence = ""
            descriptions.append(line)

    sequences.append(sequence)

    return descriptions, sequences



def parse_fasta_file(file_handle, sample_name, include_sequences=False) -> List[dict]:
    records = []

    descriptions, sequence = iterate_fasta(file_handle)
    for description in descriptions:

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

        if any(key in ERROR_COLUMNS for key in record_data.keys()):
            continue

        record_data['sample'] = sample_name
        if include_sequences:
            record_data['sequence'] = sequence

        record_data = {k: v for k,v in record_data.items() if k in COLUMNS_OF_INTEREST}

        records.append(record_data)

    return records



def read_fasta_files(base_directory, lim: Optional[int]=100):
    all_sequences = []

    count = 0
    for folder_name in os.listdir(base_directory):
        if lim is not None and count > lim:
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
    df = read_fasta_files(base_directory, lim=500)
    print(df)