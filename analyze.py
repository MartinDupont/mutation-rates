import os
import pandas as pd
import pdb
import re
from typing import List, Optional
from common import parse_fasta_file
from Levenshtein import distance
from collections import defaultdict

IDENTIFIER = 'protein'


def first_diff(string_a, string_b):
    count = 0
    for a, b in zip(string_a, string_b):
        if a != b:
            return count
        count += 1

    return None


def make_key(record):
    gene = record.get('gene')
    protein = record.get('protein')
    if not gene or not protein:
        return None

    return gene + ":" + protein


def read_and_parse_reference_genome(file_path):
    with open(file_path, 'r') as file:
        records = parse_fasta_file(file, 'REFERENCE', include_sequences=True)

    sequence_map = {}
    for r in records:
        key = make_key(r)
        sequence = r.get('sequence')
        if not key or not sequence:
            continue

        sequence_map[key] = sequence

    return sequence_map


def find_matches(records, reference_genome, mutations, anomalies, counts, anomaly_threshold=0.01):
    for r in records:
        key = make_key(r)
        sequence = r.get('sequence')
        if not key or not sequence:
            continue

        reference_sequence = reference_genome.get(key)
        if not reference_sequence:
            continue

        if sequence != reference_sequence:
            difference = distance(sequence, reference_sequence)
            mutations[key] += difference
            counts[key] += 1


            if difference > len(reference_sequence) * anomaly_threshold:
                anomalies[key] += 1




def get_diffs(base_directory, reference_genome, limit =100):

    count = 0
    mutations = defaultdict(int)
    anomalies = defaultdict(int)
    sample_counts = defaultdict(int)
    for folder_name in os.listdir(base_directory):
        if count > limit:
            break
        folder_path = os.path.join(base_directory, folder_name)
        if os.path.isdir(folder_path):
            file_path = os.path.join(folder_path, "cds_from_genomic.fna")
            if os.path.isfile(file_path):
                #print(f'reading CDS from {folder_path}')
                with open(file_path, 'r') as file:
                    parsed = parse_fasta_file(file, folder_name, include_sequences=True)
                    find_matches(parsed, reference_genome, mutations, anomalies, sample_counts)
                    count += 1

    proteins = [k for k in mutations.keys()]
    counts = [sample_counts[k] for k in mutations.keys()]
    n_mutations = [mutations[k] for k in mutations.keys()]
    ref_lengths = [len(reference_genome[k]) for k in mutations.keys()]
    n_anomalies = [anomalies[k] for k in mutations.keys()]
    sequences = [reference_genome[k] for k in mutations.keys()]

    df = pd.DataFrame({
        "identifier": proteins,
        "sample_count": counts,
        "n_mutations": n_mutations,
        "ref_length": ref_lengths,
        "n_anomalies": n_anomalies,
        "sequence": sequences,
    })


    return df


if __name__ == '__main__':
    base_directory = '/Users/martin/Documents/data/ncbi_new/ncbi_dataset/ncbi_dataset/data'
    reference_dir = '/Users/martin/Documents/data/reference_genome/ncbi_dataset/data/GCA_000005845.2/cds_from_genomic.fna'

    reference_genome = read_and_parse_reference_genome(reference_dir)
    df = get_diffs(base_directory, reference_genome)
    print("--------------------------")
    print(len(df[df['n_anomalies'] == 0]))
    print(len(df[df['n_anomalies'] < 10]))
    print(df.sort_values('n_anomalies'))