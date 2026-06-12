import unittest

from scripts.prune_redundant_consensus import (
    pairwise_similarity,
    select_nonredundant_indices,
)


class PruneRedundantConsensusTests(unittest.TestCase):
    def test_similarity_uses_union_of_non_gap_positions(self):
        self.assertEqual(pairwise_similarity("A-", "AA"), 50.0)

    def test_greedy_selection_removes_sequences_above_threshold(self):
        sequences = [
            "A" * 100,
            "A" * 99 + "C",
            "C" * 100,
        ]

        self.assertEqual(
            select_nonredundant_indices(sequences, threshold=95.0),
            [0, 2],
        )

    def test_similarity_equal_to_threshold_is_retained(self):
        sequences = [
            "A" * 100,
            "A" * 95 + "C" * 5,
        ]

        self.assertEqual(
            select_nonredundant_indices(sequences, threshold=95.0),
            [0, 1],
        )


if __name__ == "__main__":
    unittest.main()
