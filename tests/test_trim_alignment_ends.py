import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO

from scripts.trim_alignment_ends import trim_alignment_ends


class TrimAlignmentEndsTests(unittest.TestCase):
    def write_alignment(self, directory, sequences):
        path = Path(directory) / "input.fasta"
        path.write_text(
            "".join(f">{name}\n{sequence}\n" for name, sequence in sequences.items())
        )
        return path

    def read_sequences(self, path):
        with open(path) as handle:
            return {
                record.id: str(record.seq) for record in SeqIO.parse(handle, "fasta")
            }

    def test_trims_terminal_columns_below_half_acgt_coverage(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = self.write_alignment(
                temp_dir,
                {
                    "a": "--ACGTNN",
                    "b": "NNACGT--",
                    "c": "--ACGT--",
                    "d": "AAACGT--",
                },
            )
            output_path = Path(temp_dir) / "trimmed.fasta"

            stats = trim_alignment_ends(input_path, output_path, 0.5)

            self.assertEqual(
                self.read_sequences(output_path),
                {"a": "ACGT", "b": "ACGT", "c": "ACGT", "d": "ACGT"},
            )
            self.assertEqual(stats["left_trimmed"], 2)
            self.assertEqual(stats["right_trimmed"], 2)

    def test_keeps_internal_low_coverage_columns(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = self.write_alignment(
                temp_dir,
                {
                    "a": "ACGNT",
                    "b": "ACG-T",
                    "c": "ACGNT",
                    "d": "ACG-T",
                },
            )
            output_path = Path(temp_dir) / "trimmed.fasta"

            trim_alignment_ends(input_path, output_path, 0.5)

            self.assertEqual(
                self.read_sequences(output_path),
                {
                    "a": "ACGNT",
                    "b": "ACG-T",
                    "c": "ACGNT",
                    "d": "ACG-T",
                },
            )

    def test_keeps_columns_with_exactly_half_acgt_coverage(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = self.write_alignment(
                temp_dir,
                {"a": "ACGT", "b": "ACGT", "c": "NCGT", "d": "-CGT"},
            )
            output_path = Path(temp_dir) / "trimmed.fasta"

            stats = trim_alignment_ends(input_path, output_path, 0.5)

            self.assertEqual(self.read_sequences(output_path)["a"], "ACGT")
            self.assertEqual(stats["trimmed_length"], 4)

    def test_trims_terminal_n_from_single_sequence(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = self.write_alignment(temp_dir, {"a": "NNACGTNN"})
            output_path = Path(temp_dir) / "trimmed.fasta"

            trim_alignment_ends(input_path, output_path, 0.5)

            self.assertEqual(self.read_sequences(output_path), {"a": "ACGT"})


if __name__ == "__main__":
    unittest.main()
