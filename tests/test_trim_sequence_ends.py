import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO

from scripts.trim_sequence_ends import trim_sequence_ends


class TrimSequenceEndsTests(unittest.TestCase):
    def test_removes_terminal_n_and_gaps_but_keeps_internal_missing_sites(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = Path(temp_dir) / "input.fasta"
            output_path = Path(temp_dir) / "output.fasta"
            input_path.write_text(
                ">a\n--NNAC-NGTnn--\n"
                ">b\nn-A-CGT-N\n"
            )

            stats = trim_sequence_ends(input_path, output_path)

            with output_path.open() as handle:
                sequences = {
                    record.id: str(record.seq)
                    for record in SeqIO.parse(handle, "fasta")
                }
            self.assertEqual(sequences, {"a": "AC-NGT", "b": "A-CGT"})
            self.assertEqual(stats["trimmed_records"], 2)


if __name__ == "__main__":
    unittest.main()
