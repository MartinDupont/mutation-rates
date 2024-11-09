import os
import pickle
from collections import defaultdict

import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Levenshtein import distance

from common import make_key, time_it_cumulative, TIME_STORE
from common import parse_fasta_file, translate

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
def prepare_distance_data(records, sample_ids, sequence_store, genes_to_ids, current_id):
    for r in records:
        key = make_key(r)
        sample_sequence = r.get('sequence')
        if not key or not sample_sequence:
            continue

        try:
            sample_protein = translate(sample_sequence)
        except KeyError as e:
            # Invalid gene sequence
            continue
        except ValueError as e:
            continue

        id = sequence_store.get(sample_protein)
        if id is None:
            sequence_store[sample_protein] = current_id
            genes_to_ids[key].add(current_id)

            sample_ids[r['sample']].update({key: current_id})
            current_id += 1
        else:
            genes_to_ids[key].add(id)
            sample_ids[r['sample']].update({key: id})




@time_it_cumulative
def calculate_distance_matrix(sample_ids, sequence_store, genes_to_ids):
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

    reversed_sequence_store = {v:k for k, v in sequence_store.items()}

    distances = {}

    for key, sample_ids_for_key in genes_to_ids.items():
        for i in sample_ids_for_key:
            for j in sample_ids_for_key:
                s_i = reversed_sequence_store[i]
                s_j = reversed_sequence_store[j]

                difference = distance(s_i, s_j)
                distances[(i, j)] = difference
                distances[(j, i)] = difference

    dist_matrix = np.zeros((len(sample_ids), len(sample_ids)))
    i = 0
    for s1, ids_1 in sample_ids.items():
        j = 0
        for s2, ids_2 in sample_ids.items():
            differences = []
            for key, seq_id_1 in ids_1.items():
                seq_id_2 = ids_2.get(key)
                if seq_id_2 is None:
                    continue

                difference = distances[(seq_id_1, seq_id_2)]
                differences += [difference]

            total_distance = sum(differences)
            dist_matrix[i,j] = total_distance
            j += 1
        i += 1

    # df_dist_matrix = pd.DataFrame(dist_matrix, columns = sample_ids.keys(), index=sample_ids.keys())
    pickle.dump(dist_matrix, open('dist_matrix.pkl', 'wb'))
    pickle.dump(list(sample_ids.keys()), open('names.pkl', 'wb'))

    lt = extract_lower_triangular(dist_matrix)
    names = list(sample_ids.keys())

    dist_matrix = DistanceMatrix(names, lt)
    return dist_matrix




def make_phylogenetic_tree(base_directory, limit =100000):

    count = 0
    sample_ids = defaultdict(dict)
    sequence_store = defaultdict(int)
    genes_to_ids = defaultdict(set)
    current_id = 0
    for folder_name in os.listdir(base_directory):
        if count > limit:
            break
        folder_path = os.path.join(base_directory, folder_name)
        if os.path.isdir(folder_path):
            file_path = os.path.join(folder_path, "cds_from_genomic.fna")
            if os.path.isfile(file_path):
                with open(file_path, 'r') as file:
                    parsed = parse_fasta_file(file, folder_name, include_sequences=True)
                    prepare_distance_data(parsed, sample_ids, sequence_store, genes_to_ids, current_id)
                    count += 1

    pickle.dump(sample_ids, open('sample_ids.pkl', 'wb'))
    pickle.dump(sequence_store, open('sequences.pkl', 'wb'))

    print("Finished reading files")
    dist_matrix = calculate_distance_matrix(sample_ids, sequence_store, genes_to_ids)
    print(f"Finished constructing distance matrix")

    constructor = DistanceTreeConstructor()
    upgma_tree = constructor.upgma(dist_matrix)
    #NJTree = constructor.nj(distMatrix)

    Phylo.draw(upgma_tree)



if __name__ == '__main__':
    base_directory = '/Users/martin/Documents/data/ncbi_new/ncbi_dataset/ncbi_dataset/data'


    dm = make_phylogenetic_tree(base_directory)
    for key, value in TIME_STORE.items():
        print(f"\t Time for {key}: {value}")
