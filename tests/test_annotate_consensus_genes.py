#!/usr/bin/env python3

import importlib.util
import unittest
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


SCRIPT = Path(__file__).parents[1] / "scripts" / "annotate_consensus_genes.py"
SPEC = importlib.util.spec_from_file_location("annotate_genes", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class AnnotationTests(unittest.TestCase):
    def test_missing_vp2_is_inferred_from_orf3(self):
        sequence = "A" * 119 + "ATG" + "AAA" * 199 + "TAG" + "A" * 100
        reference = SeqRecord(Seq(sequence), id="GII.11_TEST")
        intervals = [
            {
                "protein": "VP1",
                "start": 20,
                "end": 100,
                "strand": "+",
                "source": "test",
            }
        ]

        inferred = MODULE.infer_vp2_interval(reference, intervals)

        self.assertEqual(inferred["start"], 120)
        self.assertEqual(inferred["end"], 722)
        self.assertEqual(inferred["source"], "inferred_from_orf3")

    def test_existing_vp2_is_not_replaced(self):
        reference = SeqRecord(Seq("A" * 1000), id="GII.11_TEST")
        intervals = [
            {
                "protein": "VP2",
                "start": 700,
                "end": 900,
                "strand": "+",
                "source": "test",
            }
        ]

        self.assertIsNone(MODULE.infer_vp2_interval(reference, intervals))

    def test_n_occupies_a_coordinate_but_gap_does_not(self):
        ref_to_cons, cons_to_ref = MODULE.build_position_map("ACN-T", "ACGT-")

        self.assertEqual(ref_to_cons, {1: 1, 2: 2, 3: 3})
        self.assertEqual(cons_to_ref, {1: 1, 2: 2, 3: 3})


if __name__ == "__main__":
    unittest.main()
