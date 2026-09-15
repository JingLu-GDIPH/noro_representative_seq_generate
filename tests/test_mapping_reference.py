#!/usr/bin/env python3

import argparse
import importlib.util
import unittest
from pathlib import Path

import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


SCRIPT = Path(__file__).parents[1] / "scripts" / "build_mapping_reference_group.py"
SPEC = importlib.util.spec_from_file_location("mapping_reference", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class MappingReferenceTests(unittest.TestCase):
    def test_distance_counts_only_acgt_in_both_sequences(self):
        self.assertAlmostEqual(MODULE.acgt_distance("ACGTN-", "ATGTAA"), 0.25)

    def test_low_posterior_site_uses_majority(self):
        args = argparse.Namespace(
            asr_pp=0.90,
            min_site_coverage=0.50,
            majority_frequency=0.60,
        )
        sequence, metrics = MODULE.hybrid_asr(
            {"a", "b", "c"},
            {"a": "AC", "b": "AT", "c": "AT"},
            np.asarray(list("AG")),
            np.asarray([0.99, 0.60]),
            args,
        )
        self.assertEqual(sequence, "AT")
        self.assertEqual(metrics["n_fraction"], 0.0)
        self.assertAlmostEqual(metrics["low_pp_fraction"], 0.5)

    def test_unresolved_low_coverage_site_becomes_n(self):
        args = argparse.Namespace(
            asr_pp=0.90,
            min_site_coverage=0.50,
            majority_frequency=0.60,
        )
        sequence, metrics = MODULE.hybrid_asr(
            {"a", "b", "c"},
            {"a": "A-", "b": "AN", "c": "AT"},
            np.asarray(list("AG")),
            np.asarray([0.99, 0.50]),
            args,
        )
        self.assertEqual(sequence, "AN")
        self.assertEqual(metrics["n_fraction"], 0.5)

    def test_snp_cluster_windows_are_collapsed(self):
        positions = [1, 2, 3, 4, 5, 6, 7, 200, 201, 202, 203, 204, 205, 206]
        self.assertEqual(MODULE.count_snp_clusters(positions), 2)

    def test_three_sequence_group_is_marked_unresolved_not_snp_rejected(self):
        records = [
            SeqRecord(Seq("ACGT"), id=f"sequence_{index}")
            for index in range(3)
        ]
        rows = MODULE.raw_qc(
            records,
            {record.id: str(record.seq) for record in records},
        )
        self.assertEqual(
            {row["flag"] for row in rows},
            {"SMALL_GROUP_NO_STABLE_REFERENCE"},
        )


if __name__ == "__main__":
    unittest.main()
