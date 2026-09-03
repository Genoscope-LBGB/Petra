import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "src" / "petra" / "workflow" / "scripts"
)
sys.path.insert(0, str(SCRIPTS_DIRECTORY))

from extract_organelles import extract_organelles

_HEADER = "sequence_id\tclass_fst_stage\tclass_snd_stage\n"
_OUTPUT_HEADER = "sequence_id\torganelle_type"


class ExtractOrganellesTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.source = self.directory / "predictions.tiara.tsv"
        self.destination = self.directory / "predictions.organelles.tsv"

    def _extract(self) -> list[str]:
        extract_organelles(self.source, self.destination)
        return self.destination.read_text(encoding="utf-8").splitlines()

    def test_lists_organelle_sequences_with_their_second_stage_type(self) -> None:
        self.source.write_text(
            _HEADER
            + "contig-1\teukarya\tn/a\n"
            + "contig-2 description\torganelle\tmitochondrion\n"
            + "contig-3\torganelle\tplastid\n"
            + "contig-4\torganelle\tunknown\n",
            encoding="utf-8",
        )

        self.assertEqual(
            self._extract(),
            [
                _OUTPUT_HEADER,
                "contig-2 description\tmitochondrion",
                "contig-3\tplastid",
                "contig-4\tunknown",
            ],
        )

    def test_writes_the_header_when_there_are_no_organelle_predictions(self) -> None:
        self.source.write_text(
            _HEADER + "contig-1\tbacteria\tn/a\n", encoding="utf-8"
        )

        self.assertEqual(self._extract(), [_OUTPUT_HEADER])

    def test_rejects_an_unexpected_tiara_table(self) -> None:
        self.source.write_text(
            "sequence_id\tclass_fst_stage\ncontig-1\torganelle\n",
            encoding="utf-8",
        )

        with self.assertRaisesRegex(ValueError, "Invalid Tiara header"):
            self._extract()


if __name__ == "__main__":
    unittest.main()
