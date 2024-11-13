import numpy as np
import pdb
import matplotlib.pyplot as plt
from Bio import Phylo
from ete3 import Tree
import os
import pickle
from common import THIS_DIRECTORY


def find_duplicates(dist_matrix, names, identity_threshold=1):
    """
    This helps us throw out samples that are clearly duplicates.
    If a duplicate pair is found, BOTH are thrown out.
    This is easier than trying to resolve clusters, like if A = B but B = C, so we need to keep only one.
    :param dist_matrix:
    :param names:
    :param identity_threshold:
    :return:
    """
    dist_matrix = np.copy(dist_matrix)
    diag = np.arange(dist_matrix.shape[0], dtype=np.uint32)
    dist_matrix[diag, diag] = np.inf

    dupes = np.argwhere(dist_matrix < identity_threshold)
    to_delete = set()
    for r in range(dupes.shape[0]):
        row = dupes[r]
        i = row[0]
        j = row[1]
        to_delete.add(names[i])
        to_delete.add(names[j])

    return list(to_delete)



if __name__ == '__main__':

    file_path = os.path.join(THIS_DIRECTORY, "dist_matrix_final_2.pkl")
    dist_matrix = pickle.load(open(file_path, 'rb'))
    file_path = os.path.join(THIS_DIRECTORY, "names.pkl")
    names = pickle.load(open(file_path, 'rb'))

    dist_matrix = dist_matrix + dist_matrix.transpose()

    to_delete = find_duplicates(dist_matrix, names)

    file_path = os.path.join(os.path.dirname(__file__), "duplicate_samples.txt")
    with open(file_path, 'w') as file:
        file.write('\n'.join(to_delete))