import pdb
import os
import pickle
from collections import defaultdict

import numpy as np
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Levenshtein import distance

from common import make_key, time_it_cumulative, TIME_STORE, read_assembly_data_report
from common import parse_fasta_file, translate, loop_over_ncbi_folders, THIS_DIRECTORY
import argparse

class IdStore():
    """
    This class just stores sequences and gives them a unique autoincrementing id.
    """
    def __init__(self):
        self.next_id = 0
        self.store = {}

    def __getitem__(self, item):
        id = self.store.get(item, None)
        if id is not None:
            return id

        retval = self.next_id
        self.store[item] = retval
        self.next_id = retval + 1
        return retval

    def to_reverse_map(self):
        return {v: k for k, v in self.store.items()}


def extract_lower_triangular(matrix):
    n = matrix.shape[0]
    lower_triangular = []

    for i in range(n):
        row = []
        for j in range(i + 1):
            row.append(matrix[i, j])
        lower_triangular.append(row)

    return lower_triangular


@time_it_cumulative
def prepare_distance_data(records, sample_ids, sequence_store: IdStore, genes_to_ids):
    for r in records:
        key = make_key(r)
        sample_id = r['sample']
        sample_sequence = r.get('sequence')
        if not key or not sample_sequence:
            continue

        if sample_ids[sample_id].get(key, None) is not None:
            # Some genomes are listed as having multiple copies of the same gene, according to NCIH
            # We filter out any extra copies.
            continue

        try:
            sample_protein = translate(sample_sequence)
        except KeyError as e:
            # Invalid gene sequence
            continue
        except ValueError as e:
            continue

        id = sequence_store[sample_protein]
        genes_to_ids[key].add(id)
        sample_ids[sample_id].update({key: id})


@time_it_cumulative
def make_reduced_distance_matrix(ids_to_sequences, genes_to_ids):
    distances = {}

    for key, sample_ids_for_key in genes_to_ids.items():
        for i in sample_ids_for_key:
            for j in sample_ids_for_key:
                if j >= i:
                    continue
                s_i = ids_to_sequences[i]
                s_j = ids_to_sequences[j]

                difference = distance(s_i, s_j)
                # This is doubling up! We should bail out if i = j
                # And put in a 0 for i == j
                distances[(i, j)] = difference
                distances[(j, i)] = difference

    return distances


@time_it_cumulative
def make_expanded_distance_matrix(sample_ids, distances):
    dist_matrix = np.zeros((len(sample_ids), len(sample_ids)))
    i = 0
    for s1, ids_1 in sample_ids.items():
        j = 0
        for s2, ids_2 in sample_ids.items():
            if j >= i:
                j += 1
                continue
            for key, seq_id_1 in ids_1.items():
                seq_id_2 = ids_2.get(key)
                if seq_id_2 is None:
                    continue

                if seq_id_1 == seq_id_2:
                    # Distance for identical pairs is 0
                    continue

                difference = distances[(seq_id_1, seq_id_2)]
                dist_matrix[i, j] += difference

            j += 1
        i += 1

    return dist_matrix


@time_it_cumulative
def calculate_distance_matrix(sample_ids, sequence_store: IdStore, genes_to_ids):
    """
    We try to calculate an NxN distance matrix.

    There are a lot of duplicated sequences. For any given gene, there are only a few variants and they are shared
    between multiple samples. So to calculate distance, it isn't efficient to stitch together the whole genome for
    each sample and compute pairwise distances.

    We instead just compute distances between unique pairs of sequences for individual genes, and for each sample and gene,
    store a reference to the sequence that it has. Due to the quadratic scaling, the runtime gains are significant.


    :param sample_ids:
    :param sequence_store:
    :param genes_to_ids:
    :return:
    """

    ids_to_sequences = sequence_store.to_reverse_map()

    distances = make_reduced_distance_matrix(ids_to_sequences, genes_to_ids)

    print(f"Done making reduced distance matrix")
    dist_matrix = make_expanded_distance_matrix(sample_ids, distances)

    return dist_matrix


def make_phylogenetic_tree(base_directory, name_to_strain, limit=100000):
    sample_ids = defaultdict(dict)
    sequence_store = IdStore()
    genes_to_ids = defaultdict(set)
    strains_encountered_so_far = set()
    for file_path, folder_name in loop_over_ncbi_folders(base_directory, limit):
        with open(file_path, 'r') as file:
            parsed = parse_fasta_file(file, folder_name, include_sequences=True)
            strain = name_to_strain.get(folder_name)
            if strain is not None and strain in strains_encountered_so_far:
                continue

            strains_encountered_so_far.add(strain)
            prepare_distance_data(parsed, sample_ids, sequence_store, genes_to_ids)

    pickle.dump(sample_ids, open(os.path.join(THIS_DIRECTORY, "sample_ids.pkl"), 'wb'))
    pickle.dump(sequence_store, open(os.path.join(THIS_DIRECTORY, "sequences.pkl"), 'wb'))

    print("Finished reading files")
    dist_matrix = calculate_distance_matrix(sample_ids, sequence_store, genes_to_ids)
    print(f"Finished constructing distance matrix")

    pickle.dump(dist_matrix, open(os.path.join(THIS_DIRECTORY, "dist_matrix.pkl"), 'wb'))
    pickle.dump(list(sample_ids.keys()), open(os.path.join(THIS_DIRECTORY, "names.pkl"), 'wb'))

    lt = extract_lower_triangular(dist_matrix)
    names = list(sample_ids.keys())

    dist_matrix = DistanceMatrix(names, lt)

    constructor = DistanceTreeConstructor()
    upgma_tree = constructor.upgma(dist_matrix)
    pickle.dump(upgma_tree, open(os.path.join(THIS_DIRECTORY, "upgma_tree.pkl"), 'wb'))

    # NJTree = constructor.nj(distMatrix)

    # Phylo.draw(upgma_tree)


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

    name_to_strain = read_assembly_data_report(args.base_directory)

    dm = make_phylogenetic_tree(args.base_directory, name_to_strain, limit=args.limit)
    print(f"Limit: {args.limit}")
    for key, value in TIME_STORE.items():
        print(f"\t Time for {key}: {value}")
