import os
import pandas as pd
import pdb
import re
from typing import List, Optional

COLUMNS_OF_INTEREST = ['gene', 'locus_tag', 'protein', 'protein_id', 'sequence']

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

    descriptions, sequences = iterate_fasta(file_handle)
    for i, description in enumerate(descriptions):

        #regex = r'\[(.*?)\]'
        regex = r'\[([^][]*(?:\[[^][]*\])*[^][]*)\]'
        matches = re.findall(regex, description)
        record_data = {}

        for match in matches:
            if '=' in match:
                key, value = match.split('=', 1)
                key = key.strip()
                value = value.strip()
                if key == 'protein_id':
                    value = value.split('.')[0]
                record_data[key] = value

        if any(key in ERROR_COLUMNS for key in record_data.keys()):
            continue

        record_data['sample'] = sample_name
        if include_sequences:
            record_data['sequence'] = sequences[i]

        record_data = {k: v for k,v in record_data.items() if k in COLUMNS_OF_INTEREST}

        records.append(record_data)

    return records