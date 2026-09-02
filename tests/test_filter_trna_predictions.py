import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "src" / "petra" / "workflow" / "scripts"
)
sys.path.insert(0, str(SCRIPTS_DIRECTORY))

from filter_trna_predictions import filter_trna_predictions  # noqa: E402

_HEADER = (
    "Sequence\t\ttRNA\tBounds\ttRNA\tAnti\tIntron Bounds\tInf\n"
    "Name    \ttRNA #\tBegin\tEnd\tType\tCodon\tBegin\tEnd\tScore\tNote\n"
    "--------\t------\t-----\t------\t----\t-----\t-----\t----\t------\t-----\n"
)
_COLUMNS = (
    "sequence_name\ttrna_number\tstart\tend\ttype\tcodon\tintron_start\t"
    "intron_end\tscore\tnote"
)


class FilterTrnaPredictionsTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.source = self.directory / "predictions.trna.raw.tsv"
        self.destination = self.directory / "predictions.trna.tsv"

    def _filter(self) -> list[str]:
        filter_trna_predictions(self.source, self.destination)
        return self.destination.read_text(encoding="utf-8").splitlines()

    def test_keeps_confident_predictions_and_names_the_columns(self) -> None:
        self.source.write_text(
            _HEADER
            + "contig-1 \t1\t500\t571\tLeu\tCAA\t0\t0\t72.4\n"
            + "contig-1 \t2\t900\t828\tGly\tGCC\t0\t0\t65.1\n",
            encoding="utf-8",
        )

        lines = self._filter()

        self.assertEqual(lines[0], _COLUMNS)
        self.assertEqual(
            [line.split("\t") for line in lines[1:]],
            [
                ["contig-1", "1", "500", "571", "Leu", "CAA", "0", "0", "72.4", "-"],
                ["contig-1", "2", "900", "828", "Gly", "GCC", "0", "0", "65.1", "-"],
            ],
        )

    def test_discards_pseudogenes_undetermined_and_suppressor_trnas(self) -> None:
        self.source.write_text(
            _HEADER
            + "contig-1\t1\t500\t571\tLeu\tCAA\t0\t0\t72.4\n"
            + "contig-1\t2\t700\t771\tPseudo\tCAA\t0\t0\t20.1\tpseudo\n"
            + "contig-1\t3\t800\t871\tUndet\t???\t0\t0\t30.2\n"
            + "contig-1\t4\t900\t971\tSup\tCTA\t0\t0\t40.3\n"
            + "contig-2\t1\t100\t171\tSeC\tTCA\t0\t0\t80.5\tIPD\n",
            encoding="utf-8",
        )

        lines = self._filter()

        self.assertEqual(
            [(line.split("\t")[0], line.split("\t")[4]) for line in lines[1:]],
            [("contig-1", "Leu"), ("contig-2", "SeC")],
        )
        self.assertEqual(lines[-1].split("\t")[-1], "IPD")

    def test_keeps_intron_coordinates_and_reversed_bounds_verbatim(self) -> None:
        self.source.write_text(
            _HEADER + "contig-1\t1\t2000\t1928\tTyr\tGTA\t1965\t1948\t60.0\n",
            encoding="utf-8",
        )

        lines = self._filter()

        self.assertEqual(
            lines[1].split("\t")[2:8], ["2000", "1928", "Tyr", "GTA", "1965", "1948"]
        )

    def test_writes_the_header_alone_without_a_prediction(self) -> None:
        self.source.write_text(_HEADER, encoding="utf-8")

        self.assertEqual(self._filter(), [_COLUMNS])

    def test_rejects_output_that_is_not_a_trnascan_table(self) -> None:
        self.source.write_text(
            "contig-1\t1\t500\t571\tLeu\tCAA\t0\t0\t72.4\n" * 4, encoding="utf-8"
        )

        with self.assertRaisesRegex(ValueError, "Unexpected tRNAscan-SE header"):
            self._filter()

    def test_rejects_a_truncated_prediction_row(self) -> None:
        self.source.write_text(
            _HEADER + "contig-1\t1\t500\t571\tLeu\n", encoding="utf-8"
        )

        with self.assertRaisesRegex(ValueError, "at least 9 fields"):
            self._filter()


if __name__ == "__main__":
    unittest.main()
