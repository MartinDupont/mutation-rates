import os
import functools
import re
import time
from collections import defaultdict
from typing import List
import json

COLUMNS_OF_INTEREST = ['gene', 'locus_tag', 'protein', 'protein_id', 'sequence', 'sample']

ERROR_COLUMNS = [
    'pseudo',
    'partial',
    'exception',
    'transl_except',
]

def make_key(record):
    gene = record.get('gene')
    protein = record.get('protein')
    if not gene or not protein:
        return None

    return gene #+ ":" + protein

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





TIME_STORE = defaultdict(lambda: 0.0)


def time_it_cumulative(func):
    """
    Decorator function to time functions
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        t = time.time()
        result = func(*args, **kwargs)

        TIME_STORE[func.__name__] += time.time() - t
        return result

    return wrapper


def time_it(func):
    """
    Decorator function to time functions
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        print(f"\nRunning {func.__name__}")
        t = time.time()
        result = func(*args, **kwargs)
        print("Time taken for %s : %4.0f sec\n" % (func.__name__, time.time() - t))
        return result

    return wrapper


TRANSLATION_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}


def translate(seq):
    protein =""
    if len(seq) % 3 != 0:
        raise ValueError("DNA sequence has length that is not a multiple of 3")

    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        protein+= TRANSLATION_TABLE[codon]
    return protein


def loop_over_ncbi_folders(base_directory, limit):
    count = 0
    for folder_name in os.listdir(base_directory):
        if count >= limit:
            break
        folder_path = os.path.join(base_directory, folder_name)
        if os.path.isdir(folder_path):
            file_path = os.path.join(folder_path, "cds_from_genomic.fna")
            if os.path.isfile(file_path):
                count += 1
                yield file_path, folder_name


def read_assembly_data_report(base_directory):
    """
    We read the assembly data because it mentions the strain name for each sample.
    We use this to filter out duplicated strains.
    :param base_directory:
    :return:
    """
    names = {}
    with open(base_directory + '/assembly_data_report.jsonl', 'r') as file:
        for line in file:
            record = json.loads(line)
            names[record['accession']] = record['organism'].get('infraspecificNames', {}).get('strain')

    return names


def get_duplicated_samples():
    file_path = os.path.join(os.path.dirname(__file__), "duplicate_samples.txt")
    if os.path.isfile(file_path):
        with open(file_path, 'r') as file:
            dupes = file.read().splitlines()
            return dupes

    raise RuntimeError("Duplicate_samples file not found")