import os
import pandas as pd
import pdb
import re
from typing import List, Optional
from common import parse_fasta_file
from Levenshtein import distance
from collections import defaultdict
from scipy.stats import poisson
import numpy as np
import matplotlib.pyplot as plt

IDENTIFIER = 'protein'


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


def translate(seq):
    protein =""
    if len(seq) % 3 != 0:
        raise ValueError("Sequence not right")

    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        protein+= TRANSLATION_TABLE[codon]
    return protein


def find_matches(records, reference_genome, mutations, anomalies, sample_counts, anomaly_threshold=0.1):
    for r in records:
        key = make_key(r)
        sample_sequence = r.get('sequence')
        if not key or not sample_sequence:
            continue

        reference_sequence = reference_genome.get(key)
        if not reference_sequence:
            continue

        reference_protein = translate(reference_sequence)
        try:
            sample_protein = translate(sample_sequence)
        except KeyError as e:
            # Invalid gene sequence
            continue

        # We count the samples in every case
        # Because we also need to track genes with no mutations to get accurate counts
        sample_counts[key] += 1

        if sample_protein != reference_protein:
            difference = distance(sample_protein, reference_protein)
            mutations[key] += difference

            if difference > len(reference_protein) * anomaly_threshold:
                anomalies[key] += 1




def get_diffs(base_directory, reference_genome, limit =100000):

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

    proteins = [k for k in sample_counts.keys()]
    counts = [sample_counts[k] for k in sample_counts.keys()]
    n_mutations = [mutations[k] for k in sample_counts.keys()]
    ref_lengths = [len(reference_genome[k]) for k in sample_counts.keys()]
    n_anomalies = [anomalies[k] for k in sample_counts.keys()]
    sequences = [reference_genome[k] for k in sample_counts.keys()]

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
    plt.ioff()
    base_directory = '/Users/martin/Documents/data/ncbi_new/ncbi_dataset/ncbi_dataset/data'
    reference_dir = '/Users/martin/Documents/data/reference_genome/ncbi_dataset/data/GCA_000005845.2/cds_from_genomic.fna'

    reference_genome = read_and_parse_reference_genome(reference_dir)
    df = get_diffs(base_directory, reference_genome)

    import pickle
    pickle.dump(df, open('df', 'wb'))
