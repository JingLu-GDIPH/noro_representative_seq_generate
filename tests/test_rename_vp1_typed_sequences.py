import importlib.util
from pathlib import Path
import unittest


SCRIPT = Path(__file__).parents[1] / "scripts" / "rename_vp1_typed_sequences.py"


def load_module():
    spec = importlib.util.spec_from_file_location("rename_vp1", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class RenameVp1TypedSequencesTests(unittest.TestCase):
    def test_safe_collection_date_normalizes_separators(self):
        module = load_module()
        self.assertEqual(
            module.safe_collection_date("12 May/2025"),
            "12-May-2025",
        )

    def test_build_seqname_uses_vp1_genotype(self):
        module = load_module()
        row = {
            "group": "GII",
            "accession": "PZ478197",
            "collection_date": "12-May-2025",
            "vp1_genotype": "GII.4",
        }
        self.assertEqual(
            module.build_seqname(row),
            "GII.4_PZ478197_12-May-2025",
        )

    def test_build_seqname_marks_no_hit_as_unassigned(self):
        module = load_module()
        row = {
            "group": "GI",
            "accession": "LC726063",
            "collection_date": "2022-01",
            "vp1_genotype": "",
        }
        self.assertEqual(
            module.build_seqname(row),
            "GI.Unassigned_LC726063_2022-01",
        )

    def test_accession_is_recovered_from_renamed_fasta_id(self):
        module = load_module()
        self.assertEqual(
            module.accession_from_fasta_id(
                "GII.4_PZ478197_12-May-2025"
            ),
            "PZ478197",
        )


if __name__ == "__main__":
    unittest.main()
