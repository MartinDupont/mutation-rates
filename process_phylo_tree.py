import os
import pickle
from collections import defaultdict

import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Levenshtein import distance

from common import make_key, time_it_cumulative, TIME_STORE, time_it
from common import parse_fasta_file, translate



class IdStore():
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
        return {v:k for k, v in self.store.items()}



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
def make_reduced_distance_matrix(reversed_sequence_store, genes_to_ids):

    distances = {}

    for key, sample_ids_for_key in genes_to_ids.items():
        for i in sample_ids_for_key:
            for j in sample_ids_for_key:
                if j >= i:
                    continue
                s_i = reversed_sequence_store[i]
                s_j = reversed_sequence_store[j]

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
            for key, seq_id_1 in ids_1.items():
                seq_id_2 = ids_2.get(key)
                if seq_id_2 is None:
                    continue

                if seq_id_1 == seq_id_2:
                    # Distance for identical pairs is 0
                    continue

                difference = distances[(seq_id_1, seq_id_2)]
                dist_matrix[i,j] += difference

            j += 1
        i += 1

    return dist_matrix


@time_it_cumulative
def calculate_distance_matrix(sample_ids, sequence_store:IdStore, genes_to_ids):
    """
    We try to calculate an NxN distance matrix.

    There are a lot of duplicated sequences. For any given gene, there are only a few variants and they are shared
    between multiple samples. So to calculate distance, it isn't efficient to stitch together the whole genome for
    each sample and compute pairwise distances.

    We instead just compute distances between pairs of sequences for individual genes, and for each sample and gene,
    store a reference to the sequence that it has. Due to the quadratic scaling, the runtime gains are significant.


    :param sample_ids:
    :param sequence_store:
    :param genes_to_ids:
    :return:
    """

    reversed_sequence_store = sequence_store.to_reverse_map()

    distances = make_reduced_distance_matrix(reversed_sequence_store, genes_to_ids)

    dist_matrix = make_expanded_distance_matrix(sample_ids, distances)

    return dist_matrix




def make_phylogenetic_tree(base_directory, limit =100000):

    count = 0
    sample_ids = defaultdict(dict)
    sequence_store = IdStore()
    genes_to_ids = defaultdict(set)
    for folder_name in os.listdir(base_directory):
        if count >= limit:
            break
        folder_path = os.path.join(base_directory, folder_name)
        if os.path.isdir(folder_path):
            file_path = os.path.join(folder_path, "cds_from_genomic.fna")
            if os.path.isfile(file_path):
                with open(file_path, 'r') as file:
                    parsed = parse_fasta_file(file, folder_name, include_sequences=True)
                    prepare_distance_data(parsed, sample_ids, sequence_store, genes_to_ids)
                    count += 1

    pickle.dump(sample_ids, open('sample_ids.pkl', 'wb'))
    pickle.dump(sequence_store, open('sequences.pkl', 'wb'))

    print("Finished reading files")
    dist_matrix = calculate_distance_matrix(sample_ids, sequence_store, genes_to_ids)
    print(f"Finished constructing distance matrix")

    pickle.dump(dist_matrix, open('dist_matrix.pkl', 'wb'))
    pickle.dump(list(sample_ids.keys()), open('names.pkl', 'wb'))

    lt = extract_lower_triangular(dist_matrix)
    names = list(sample_ids.keys())

    dist_matrix = DistanceMatrix(names, lt)

    constructor = DistanceTreeConstructor()
    upgma_tree = constructor.upgma(dist_matrix)
    pickle.dump(upgma_tree, open('upgma_tree.pkl', 'wb'))

    #NJTree = constructor.nj(distMatrix)

    #Phylo.draw(upgma_tree)



if __name__ == '__main__':
    base_directory = '/Users/martin/Documents/data/ncbi_new/ncbi_dataset/ncbi_dataset/data'

    limit = 1000000
    TIME_STORE.clear()
    dm = make_phylogenetic_tree(base_directory, limit=limit)
    print(f"Limit: {limit}")
    for key, value in TIME_STORE.items():
        print(f"\t Time for {key}: {value}")
