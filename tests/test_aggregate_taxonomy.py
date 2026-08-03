import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "src" / "petra" / "workflow" / "scripts"
)
sys.path.insert(0, str(SCRIPTS_DIRECTORY))

from aggregate_taxonomy import aggregate_taxonomy  # noqa: E402
from assign_taxonomy import TAXONOMIC_RANKS  # noqa: E402


class AggregateTaxonomyTest(unittest.TestCase):
    def test_combines_metadata_reasons_and_whokaryote_disagreements(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            raw_taxonomy = directory / "raw.tsv"
            coverage = directory / "coverage.tsv"
            assignment_reasons = directory / "assignment_reasons.tsv"
            whokaryote_reasons = directory / "whokaryote_reasons.tsv"
            destination = directory / "final.tsv"
            self._write_taxonomy(
                raw_taxonomy,
                [
                    self._taxonomy_row("assigned", "Bacteria"),
                    self._taxonomy_row("low_density", "N/A"),
                    self._taxonomy_row("disagreement", "Bacteria"),
                ],
            )
            self._write_coverages(
                coverage,
                [
                    ("disagreement", "3"),
                    ("assigned", "10.5"),
                    ("low_density", "0"),
                ],
            )
            self._write_reasons(
                assignment_reasons,
                [
                    ("assigned", "-"),
                    ("low_density", "gene_density"),
                    ("disagreement", "-"),
                ],
            )
            self._write_reasons(
                whokaryote_reasons,
                [
                    ("assigned", "-"),
                    ("low_density", "-"),
                    ("disagreement", "whokaryote"),
                ],
            )

            aggregate_taxonomy(
                raw_taxonomy,
                coverage,
                assignment_reasons,
                whokaryote_reasons,
                destination,
            )

            rows = [line.split("\t") for line in destination.read_text().splitlines()]
            self.assertEqual(
                rows[0],
                [
                    "sequence_name",
                    "sequence_length",
                    "median_coverage",
                    "number_of_aligned_clusters",
                    "number_of_genes_in_aligned_clusters",
                    *TAXONOMIC_RANKS,
                    "unassigned_reason",
                ],
            )
            self.assertEqual(
                rows[1][0:6], ["assigned", "100", "10.5", "1", "1", "Bacteria"]
            )
            self.assertEqual(rows[1][-1], "-")
            self.assertEqual(
                rows[2][0:6], ["low_density", "100", "0", "1", "1", "N/A"]
            )
            self.assertEqual(rows[2][-1], "gene_density")
            self.assertEqual(
                rows[3],
                [
                    "disagreement",
                    "100",
                    "3",
                    "1",
                    "1",
                    *(["N/A"] * len(TAXONOMIC_RANKS)),
                    "whokaryote",
                ],
            )

    def test_rejects_multiple_reasons_for_one_sequence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            raw_taxonomy = directory / "raw.tsv"
            coverage = directory / "coverage.tsv"
            assignment_reasons = directory / "assignment_reasons.tsv"
            whokaryote_reasons = directory / "whokaryote_reasons.tsv"
            self._write_taxonomy(
                raw_taxonomy, [self._taxonomy_row("contig", "Bacteria")]
            )
            self._write_coverages(coverage, [("contig", "5")])
            self._write_reasons(assignment_reasons, [("contig", "no_consensus")])
            self._write_reasons(whokaryote_reasons, [("contig", "whokaryote")])

            with self.assertRaisesRegex(ValueError, "conflicting unassigned reasons"):
                aggregate_taxonomy(
                    raw_taxonomy,
                    coverage,
                    assignment_reasons,
                    whokaryote_reasons,
                    directory / "final.tsv",
                )

    @staticmethod
    def _taxonomy_row(sequence_id: str, domain: str) -> list[str]:
        return [
            sequence_id,
            "100",
            "1",
            "1",
            domain,
            *(["N/A"] * (len(TAXONOMIC_RANKS) - 1)),
        ]

    @staticmethod
    def _write_taxonomy(destination: Path, rows: list[list[str]]) -> None:
        header = [
            "sequence_name",
            "sequence_length",
            "number_of_aligned_clusters",
            "number_of_genes_in_aligned_clusters",
            *TAXONOMIC_RANKS,
        ]
        destination.write_text(
            "\n".join("\t".join(row) for row in [header, *rows]) + "\n",
            encoding="utf-8",
        )

    @staticmethod
    def _write_coverages(destination: Path, rows: list[tuple[str, str]]) -> None:
        destination.write_text(
            "sequence_name\tmedian_coverage\n"
            + "".join(
                f"{sequence_id}\t{median_coverage}\n"
                for sequence_id, median_coverage in rows
            ),
            encoding="utf-8",
        )

    @staticmethod
    def _write_reasons(destination: Path, rows: list[tuple[str, str]]) -> None:
        destination.write_text(
            "sequence_name\tunassigned_reason\n"
            + "".join(f"{sequence_id}\t{reason}\n" for sequence_id, reason in rows),
            encoding="utf-8",
        )


if __name__ == "__main__":
    unittest.main()
