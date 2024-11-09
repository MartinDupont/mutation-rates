import numpy as np
from numpy.testing import assert_array_equal

from process_phylo_tree import make_reduced_distance_matrix, make_expanded_distance_matrix


def test_make_reduced_distance_matrix():
    # There are three genes, one has three variants, one has two variants, one has no variants
    #
    genes_to_ids = {
        'gene_1': {0, 1, 2},
        'gene_2': {3, 4},
        'gene_3': {5}
    }

    ids_to_sequences = {
        0: 'AAAAA',
        1: 'AAAAB',
        2: 'AAABB',
        3: 'CCCCC',
        4: 'DDDDD'
    }

    expected_result = {
        # gene 1
        (0, 1): 1, # Sequence 0 and 1 differ at a single location
        (1, 0): 1,
        (0, 2): 2, # Sequence 0 and 2 differ at two locations
        (2, 0): 2,
        (1, 2): 1, # Sequence 1 and 2 differ at a single location
        (2, 1): 1,
        # gene 2
        (3, 4): 5, # Sequence 2 and 3 differ at every location
        (4, 3): 5
        # we expect no entries for gene 3 as there are no pairs
    }

    result = make_reduced_distance_matrix(ids_to_sequences, genes_to_ids)

    assert expected_result == result



def test_make_expanded_distance_matrix():
    # Same setup as above
    # There are three genes, one has three variants, one has two variants, one has no variants

    reduced_distance_matrix = {
        # gene 1. For convenience we say that all variants are equally far apart
        (0, 1): 1,
        (1, 0): 1,
        (0, 2): 1,
        (2, 0): 1,
        (1, 2): 1,
        (2, 1): 1,
        # gene 2, they differ by two
        (3, 4): 2,
        (4, 3): 2,
        # gene 3 has no variants

        # Irrelevant combinations should be ignored
        (42, 103): 42
    }

    sample_ids = {
        'sample_1': { 'gene_1': 0, 'gene_2': 3, 'gene_3': 5},
        'sample_2': { 'gene_1': 1, 'gene_2': 4, 'gene_3': 5},
        'sample_3': { 'gene_1': 2, 'gene_2': 3}, # missing genes should be handled
    }


    expected_result = np.array(
        [
            # sample 1 differs from sample 2 at both genes, so the difference is 1 + 2 = 3
            # meanwhile sample 1 differs from sample 3 at only gene 1, so the distance is 1
            [0, 3, 1],
            # sample 2 differs from sample 3 at gene_1 and gene_2
            [3, 0, 3],
            [1, 3, 0],
        ]
    )

    result = make_expanded_distance_matrix(sample_ids, reduced_distance_matrix)

    assert_array_equal(result, expected_result)
