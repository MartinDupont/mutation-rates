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


def find_matches(records, sample_counts, unique_sequences, sequence_counts, reference_genome):
    for r in records:
        key = make_key(r)
        sample_sequence = r.get('sequence')
        if not key or not sample_sequence:
            continue

        # We don't actually need to compare to the reference sequence,
        # But we only calculate mutations for genes that are found on the reference sequence.
        reference_sequence = reference_genome.get(key)
        if not reference_sequence:
            continue

        try:
            sample_protein = translate(sample_sequence)
        except KeyError as e:
            # Invalid gene sequence
            continue

        sample_counts[key] += 1
        unique_sequences[key].add(sample_protein)
        sequence_counts[sample_protein] += 1


def count_mutations_maf(sample_counts, unique_sequences, sequence_counts, anomalies, anomaly_threshold=0.1):
    """
    The reference genome isn't necessarily the one with the major allele frequency.
    So, we store all the unique sequences encountered and how many times each sequence has been encountered.
    We take the most common sequence and calculate the number of mutations by summing the edit distance
    between all encountered sequences and the most common one, weighted by  occurrence.

    This is equivalent to summing mutations on a per-allele basis if the most common sequence is >50%
    Empirically, this is true > 90% of the time. Good enough for now.

    :param sample_counts:
    :param unique_sequences:
    :param sequence_counts:
    :return:
    """

    mutations_maf = defaultdict(int)
    mafs = {}
    for key, value in sample_counts.items():
        sequences = unique_sequences[key]
        biggest_count = 0
        most_common_sequence = None
        for s in sequences:
            if sequence_counts[s] > biggest_count:
                biggest_count = sequence_counts[s]
                most_common_sequence = s

        for s in sequences:
            difference = distance(most_common_sequence, s)
            mutations_maf[key] += difference * sequence_counts[s]

            if difference > len(most_common_sequence) * anomaly_threshold:
                anomalies[key] += 1

        mafs[key] = sequence_counts[most_common_sequence] /  sum(sequence_counts[s] for s in sequences)


    return mutations_maf, mafs



def get_diffs(base_directory, reference_genome, limit =100000):

    count = 0
    anomalies = defaultdict(int)
    sample_counts = defaultdict(int)
    unique_sequences = defaultdict(set)
    sequence_counts = defaultdict(int)
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
                    find_matches(parsed, sample_counts, unique_sequences, sequence_counts, reference_genome)
                    count += 1


    mutations, mafs = count_mutations_maf(sample_counts, unique_sequences, sequence_counts, anomalies)

    proteins = [k for k in sample_counts.keys()]
    counts = [sample_counts[k] for k in sample_counts.keys()]
    n_mutations = [mutations[k] for k in sample_counts.keys()]
    ref_lengths = [len(reference_genome[k]) for k in sample_counts.keys()]
    n_anomalies = [anomalies[k] for k in sample_counts.keys()]
    n_sequences = [len(unique_sequences[k]) for k in sample_counts.keys()]
    #sequences = [reference_genome[k] for k in sample_counts.keys()]
    major_alle_freqs = [mafs[k] for k in sample_counts.keys()]

    df = pd.DataFrame({
        "identifier": proteins,
        "sample_count": counts,
        "n_mutations": n_mutations,
        "ref_length": ref_lengths,
        "n_anomalies": n_anomalies,
        "n_sequences": n_sequences,
        "major_allele_freq": major_alle_freqs,
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
