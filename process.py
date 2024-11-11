import os
import pdb
from collections import defaultdict

import matplotlib.pyplot as plt
import pandas as pd
from Levenshtein import distance

from common import parse_fasta_file, make_key, translate, loop_over_ncbi_folders


def first_diff(string_a, string_b):
    count = 0
    for a, b in zip(string_a, string_b):
        if a != b:
            return count
        count += 1

    return None



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

def find_matches(records, sample_counts, unique_sequences, sequence_counts, reference_genome, protein_names):
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
        except ValueError as e:
            # Invalid gene sequence
            continue

        sample_counts[key] += 1
        unique_sequences[key].add(sample_protein)
        sequence_counts[sample_protein] += 1
        # This will clobber any existing key. We don't really care because the protein names are just for flavor
        protein_names[key] = r.get('protein')


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

    anomalies = defaultdict(int)
    sample_counts = defaultdict(int)
    unique_sequences = defaultdict(set)
    sequence_counts = defaultdict(int)
    protein_names = dict()
    for file_path, folder_name in loop_over_ncbi_folders(base_directory, limit):
        with open(file_path, 'r') as file:
            parsed = parse_fasta_file(file, folder_name, include_sequences=True)
            find_matches(parsed, sample_counts, unique_sequences, sequence_counts, reference_genome, protein_names)


    mutations, mafs = count_mutations_maf(sample_counts, unique_sequences, sequence_counts, anomalies)

    identifiers = [k for k in sample_counts.keys()]
    proteins = [protein_names[k] for k in sample_counts.keys()]
    counts = [sample_counts[k] for k in sample_counts.keys()]
    n_mutations = [mutations[k] for k in sample_counts.keys()]
    ref_lengths = [len(reference_genome[k]) for k in sample_counts.keys()]
    n_anomalies = [anomalies[k] for k in sample_counts.keys()]
    n_sequences = [len(unique_sequences[k]) for k in sample_counts.keys()]
    #sequences = [reference_genome[k] for k in sample_counts.keys()]
    major_alle_freqs = [mafs[k] for k in sample_counts.keys()]

    df = pd.DataFrame({
        "identifier": identifiers,
        "protein": proteins,
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
