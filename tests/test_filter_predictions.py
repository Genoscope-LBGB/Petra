import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "src" / "petra" / "workflow" / "scripts"
)
sys.path.insert(0, str(SCRIPTS_DIRECTORY))

from filter_predictions import filter_predictions  # noqa: E402

TABLE = (
    "sequence_name\tsequence_length\tgc_percent\tmedian_coverage\tphylum\tspecies\n"
    "contig-a\t100\t60.5\t20\tActinopteri\tChelon auratus\n"
    "contig-b\t200\t61.0\t50\tChordata\tN/A\n"
    "contig-c\t300\t58.0\t30\tProteobacteria\tN/A\n"
    "contig-d\t400\t62.0\tN/A\tN/A\tN/A\n"
)


class FilterPredictionsTest(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self._temporary_directory.cleanup)
        directory = Path(self._temporary_directory.name)
        self.source = directory / "predictions.final.tsv"
        self.source.write_text(TABLE, encoding="utf-8")
        self.destination = directory / "kept.tsv"

    def _retained(self, expression: str) -> list[str]:
        filter_predictions(self.source, self.destination, expression)
        lines = self.destination.read_text(encoding="utf-8").splitlines()
        self.assertEqual(lines[0].split("\t")[0], "sequence_name")
        return [line.split("\t")[0] for line in lines[1:]]

    def test_combines_conditions_with_and_and_not(self) -> None:
        self.assertEqual(
            self._retained(
                "median_coverage<40 and gc_percent>59 and not phylum==Actinopteri"
            ),
            [],
        )
        self.assertEqual(
            self._retained("median_coverage<40 and gc_percent>59"), ["contig-a"]
        )

    def test_reads_not_written_as_an_exclamation_mark(self) -> None:
        self.assertEqual(
            self._retained("!phylum==Actinopteri and median_coverage<40"), ["contig-c"]
        )

    def test_or_binds_looser_than_and(self) -> None:
        self.assertEqual(
            self._retained("phylum==Chordata or sequence_length<400 and gc_percent<59"),
            ["contig-b", "contig-c"],
        )

    def test_parentheses_regroup_conditions(self) -> None:
        self.assertEqual(
            self._retained(
                "(phylum==Chordata or sequence_length<400) and gc_percent<59"
            ),
            ["contig-c"],
        )

    def test_quotes_hold_the_spaces_of_a_value(self) -> None:
        self.assertEqual(self._retained('species=="Chelon auratus"'), ["contig-a"])
        self.assertEqual(self._retained("species=='Chelon auratus'"), ["contig-a"])

    def test_keywords_are_case_insensitive(self) -> None:
        self.assertEqual(
            self._retained("gc_percent>59 AND NOT phylum==Actinopteri"),
            ["contig-b", "contig-d"],
        )

    def test_non_numeric_value_never_satisfies_a_numeric_comparison(self) -> None:
        self.assertEqual(
            self._retained("median_coverage>=0"), ["contig-a", "contig-b", "contig-c"]
        )
        self.assertEqual(
            self._retained("median_coverage!=20"),
            ["contig-b", "contig-c", "contig-d"],
        )

    def test_matches_text_exactly(self) -> None:
        self.assertEqual(self._retained("phylum==N/A"), ["contig-d"])
        self.assertEqual(self._retained("phylum==actinopteri"), [])

    def test_rejects_an_attribute_the_table_lacks(self) -> None:
        with self.assertRaisesRegex(ValueError, "Unknown attribute.*gc_percentage"):
            filter_predictions(self.source, self.destination, "gc_percentage>59")

    def test_rejects_ordering_against_text(self) -> None:
        with self.assertRaisesRegex(ValueError, "non-numeric value 'Chordata'"):
            filter_predictions(self.source, self.destination, "phylum>Chordata")

    def test_rejects_a_word_that_is_not_a_condition(self) -> None:
        with self.assertRaisesRegex(ValueError, "'median_coverage' is not a condition"):
            filter_predictions(self.source, self.destination, "median_coverage")

    def test_rejects_an_unbalanced_expression(self) -> None:
        with self.assertRaisesRegex(ValueError, "never closed by a"):
            filter_predictions(self.source, self.destination, "(gc_percent>59")
        with self.assertRaisesRegex(ValueError, "no matching"):
            filter_predictions(self.source, self.destination, "gc_percent>59)")
        with self.assertRaisesRegex(ValueError, "condition is missing"):
            filter_predictions(self.source, self.destination, "gc_percent>59 and")
        with self.assertRaisesRegex(ValueError, "no condition on its left"):
            filter_predictions(self.source, self.destination, "or gc_percent>59")
        with self.assertRaisesRegex(ValueError, "after a complete expression"):
            filter_predictions(self.source, self.destination, "gc_percent>59 phylum==N/A")

    def test_rejects_an_unclosed_quote(self) -> None:
        with self.assertRaisesRegex(ValueError, "quote is never closed"):
            filter_predictions(self.source, self.destination, "species=='Chelon")

    def test_rejects_an_empty_expression(self) -> None:
        with self.assertRaisesRegex(ValueError, "holds no condition"):
            filter_predictions(self.source, self.destination, "   ")

    def test_rejects_a_row_with_a_wrong_column_count(self) -> None:
        self.source.write_text(TABLE + "contig-e\t500\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "line 6: expected 6 columns, found 2"):
            filter_predictions(self.source, self.destination, "gc_percent>0")


if __name__ == "__main__":
    unittest.main()
