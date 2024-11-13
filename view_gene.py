import os
from collections import defaultdict
import argparse
import matplotlib.pyplot as plt
from Levenshtein import distance

from common import parse_fasta_file, make_key, translate, loop_over_ncbi_folders


def make_muations_list(unique_sequences, sequence_counts):
    biggest_count = 0
    most_common_sequence = None
    for s in unique_sequences:
        if sequence_counts[s] > biggest_count:
            biggest_count = sequence_counts[s]
            most_common_sequence = s

    n_mutations_list = []

    for s in unique_sequences:
        difference = distance(most_common_sequence, s)
        for i in range(sequence_counts[s]):
            n_mutations_list += [difference]

    return n_mutations_list


def find_matches(target_key, records, unique_sequences, sequence_counts):
    for r in records:
        key = make_key(r)
        sample_sequence = r.get('sequence')
        if not key or not sample_sequence:
            continue

        if key != target_key:
            continue

        try:
            sample_protein = translate(sample_sequence)
        except KeyError as e:
            # Invalid gene sequence
            continue

        unique_sequences.add(sample_protein)
        sequence_counts[sample_protein] += 1
        break


def view_gene(base_directory, target_key, limit=100000):
    unique_sequences = set()
    sequence_counts = defaultdict(int)
    for file_path, folder_name in loop_over_ncbi_folders(base_directory, limit):
        with open(file_path, 'r') as file:
            parsed = parse_fasta_file(file, folder_name, include_sequences=True)
            find_matches(target_key, parsed, unique_sequences, sequence_counts)

    n_mutations_list = make_muations_list(unique_sequences, sequence_counts)

    plt.hist(n_mutations_list, 20, alpha=0.5, color="blue")
    plt.xlabel('N differences from most common sequence')
    plt.ylabel('Count')
    plt.title(target_key)
    plt.show()
    plt.close()


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
        "--gene",
        "-g",
        dest="key",
        help="The ID of the gene",
        type=str,
    )

    args = parser.parse_args()
    print(args)

    plt.ioff()

    key = 'rhsC:rhs element protein RhsC'
    key = 'dgoK:2-dehydro-3-deoxygalactonokinase'
    key = 'hofN:DNA utilization protein HofN'
    key = 'nanK:N-acetylmannosamine kinase'

    view_gene(args.base_directory, args.key)
