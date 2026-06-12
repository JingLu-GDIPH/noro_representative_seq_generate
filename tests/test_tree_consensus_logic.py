import tempfile
import unittest
from pathlib import Path

from Bio import Phylo
from io import StringIO

from scripts import enhanced_consensus_with_validation as consensus
from scripts import generate_consensus


class TreeConsensusLogicTests(unittest.TestCase):
    def test_enhanced_genotype_extraction_accepts_gi_from_filename(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            fasta = Path(temp_dir) / "GI.3_consensus_2.fasta"
            fasta.write_text(">node_1\nACGT\n")

            self.assertEqual(
                consensus.extract_common_genotype_from_cluster(fasta),
                "GI.3",
            )

    def test_enhanced_genotype_extraction_accepts_gi_from_sequence_ids(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            fasta = Path(temp_dir) / "iteration_input.fasta"
            fasta.write_text(">GI.6_accession_date\nACGT\n")

            self.assertEqual(
                consensus.extract_common_genotype_from_cluster(fasta),
                "GI.6",
            )

    def test_enhanced_node_genotypes_accept_gi(self):
        tree = Phylo.read(
            StringIO("(GI.7_accession_a:1,GI.7_accession_b:1)Root;"),
            "newick",
        )

        self.assertEqual(
            set(
                consensus.get_node_genotypes(
                    tree,
                    tree.root,
                    {},
                    {
                        "GI.7_accession_a": "AAAA",
                        "GI.7_accession_b": "AAAT",
                    },
                )
            ),
            {"GI.7"},
        )

    def test_parse_iqtree_state_reconstructs_node_sequences(self):
        state_text = """\
# IQ-TREE ancestral states
Node\tSite\tState\tp_A\tp_C\tp_G\tp_T
Node1\t1\tA\t0.9\t0.1\t0.0\t0.0
Node2\t1\tG\t0.0\t0.0\t0.9\t0.1
Node1\t2\tC\t0.1\t0.9\t0.0\t0.0
Node2\t2\tT\t0.0\t0.0\t0.1\t0.9
"""
        with tempfile.TemporaryDirectory() as temp_dir:
            state_file = Path(temp_dir) / "tree.state"
            state_file.write_text(state_text)
            self.assertEqual(
                consensus.parse_ancestral_sequences(state_file),
                {"Node1": "AC", "Node2": "GT"},
            )

    def test_clade_sequences_exclude_ancestral_sequence(self):
        tree = Phylo.read(StringIO("((a:1,b:1)N1:1,c:1)Root;"), "newick")
        node = next(clade for clade in tree.get_nonterminals() if clade.name == "N1")
        sequences = consensus.get_node_sequences(
            tree,
            node,
            {"N1": "TTTT"},
            {"a": "AAAA", "b": "AAAT", "c": "CCCC"},
        )
        self.assertEqual(sequences, ["AAAA", "AAAT"])

    def test_optimal_nodes_form_complete_non_overlapping_partition(self):
        tree = Phylo.read(StringIO("((a:1,b:1)N1:1,c:1)Root;"), "newick")
        selected = generate_consensus.find_optimal_nodes(
            tree,
            {},
            {"a": "A" * 100, "b": "A" * 96 + "C" * 4, "c": "C" * 100},
            95.0,
        )
        selected_tip_sets = [
            {tip.name for tip in node.get_terminals()}
            for node, _, _ in selected
        ]
        self.assertEqual(selected_tip_sets, [{"a", "b"}, {"c"}])
        flattened = [tip for tip_set in selected_tip_sets for tip in tip_set]
        self.assertEqual(sorted(flattened), ["a", "b", "c"])

    def test_root_is_selected_only_once(self):
        tree = Phylo.read(StringIO("(a:1,b:1)Root;"), "newick")
        selected = generate_consensus.find_optimal_nodes(
            tree,
            {},
            {"a": "A" * 100, "b": "A" * 96 + "C" * 4},
            95.0,
        )
        self.assertEqual(len(selected), 1)
        self.assertIs(selected[0][0], tree.root)


if __name__ == "__main__":
    unittest.main()
