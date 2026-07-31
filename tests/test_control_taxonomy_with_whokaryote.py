import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "src" / "petra" / "workflow" / "scripts"
)
sys.path.insert(0, str(SCRIPTS_DIRECTORY))

from assign_taxonomy import TAXONOMIC_RANKS  # noqa: E402
from control_taxonomy_with_whokaryote import control_taxonomy  # noqa: E402


class ControlTaxonomyWithWhokaryoteTest(unittest.TestCase):
    def test_nulls_disagreements_and_preserves_other_assignments(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            raw_taxonomy = directory / "raw.tsv"
            predictions = directory / "whokaryote.tsv"
            destination = directory / "controlled.tsv"
            raw_rows = [
                self._raw_row("eukaryote", "3", "Eukaryota", kingdom="Fungi"),
                self._raw_row("bacterium", "2", "Bacteria"),
                self._raw_row("archaeon", "4", "Archaea"),
                self._raw_row("virus", "1", "Viruses"),
                self._raw_row("viroid", "1", "Viroids"),
                self._raw_row("not_predicted", "2", "Bacteria"),
                self._raw_row(
                    "petra_na", "0", "N/A", reason="gene_density", kingdom="Fungi"
                ),
            ]
            self._write_raw(raw_taxonomy, raw_rows)
            predictions.write_text(
                "contig\tpredicted\n"
                "eukaryote\teukaryote\n"
                "bacterium\teukaryote\n"
                "archaeon\tprokaryote\n"
                "virus\tprokaryote\n"
                "viroid\teukaryote\n"
                "petra_na\teukaryote\n",
                encoding="utf-8",
            )

            control_taxonomy(raw_taxonomy, predictions, destination)

            rows = [line.split("\t") for line in destination.read_text().splitlines()]
            self.assertEqual(
                rows[0],
                [
                    "sequence_name",
                    "number_of_aligned_clusters",
                    "number_of_genes_in_aligned_clusters",
                    *TAXONOMIC_RANKS,
                    "unassigned_reason",
                ],
            )
            self.assertEqual(rows[1], raw_rows[0])
            self.assertEqual(
                rows[2],
                [
                    "bacterium",
                    "2",
                    "2",
                    *(["N/A"] * len(TAXONOMIC_RANKS)),
                    "whokaryote",
                ],
            )
            self.assertEqual(rows[3], raw_rows[2])
            self.assertEqual(
                rows[4],
                [
                    "virus",
                    "1",
                    "1",
                    *(["N/A"] * len(TAXONOMIC_RANKS)),
                    "whokaryote",
                ],
            )
            self.assertEqual(
                rows[5],
                [
                    "viroid",
                    "1",
                    "1",
                    *(["N/A"] * len(TAXONOMIC_RANKS)),
                    "whokaryote",
                ],
            )
            self.assertEqual(rows[6], raw_rows[5])
            self.assertEqual(rows[7], raw_rows[6])

    def test_rejects_unknown_whokaryote_prediction(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            raw_taxonomy = directory / "raw.tsv"
            predictions = directory / "whokaryote.tsv"
            self._write_raw(raw_taxonomy, [self._raw_row("contig", "1", "Bacteria")])
            predictions.write_text(
                "contig\tpredicted\ncontig\tvirus\n", encoding="utf-8"
            )

            with self.assertRaisesRegex(ValueError, "Unexpected Whokaryote prediction"):
                control_taxonomy(
                    raw_taxonomy, predictions, directory / "controlled.tsv"
                )

    def test_rejects_unknown_whokaryote_identifier(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            raw_taxonomy = directory / "raw.tsv"
            predictions = directory / "whokaryote.tsv"
            self._write_raw(raw_taxonomy, [self._raw_row("contig", "1", "Bacteria")])
            predictions.write_text(
                "contig\tpredicted\nother\tprokaryote\n", encoding="utf-8"
            )

            with self.assertRaisesRegex(ValueError, "absent from the raw taxonomy"):
                control_taxonomy(
                    raw_taxonomy, predictions, directory / "controlled.tsv"
                )

    def test_rejects_duplicate_identifiers(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            raw_taxonomy = directory / "raw.tsv"
            predictions = directory / "whokaryote.tsv"
            self._write_raw(raw_taxonomy, [self._raw_row("contig", "1", "Bacteria")])
            predictions.write_text(
                "contig\tpredicted\n" "contig\tprokaryote\n" "contig\teukaryote\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "Duplicate Whokaryote"):
                control_taxonomy(
                    raw_taxonomy, predictions, directory / "controlled.tsv"
                )

    def test_rejects_unknown_petra_domain(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            raw_taxonomy = directory / "raw.tsv"
            predictions = directory / "whokaryote.tsv"
            self._write_raw(raw_taxonomy, [self._raw_row("contig", "1", "Other")])
            predictions.write_text("contig\tpredicted\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "Unexpected PETRA domain"):
                control_taxonomy(
                    raw_taxonomy, predictions, directory / "controlled.tsv"
                )

    @staticmethod
    def _raw_row(
        sequence_id: str,
        aligned_clusters: str,
        domain: str,
        reason: str = "-",
        **ranks: str,
    ) -> list[str]:
        lineage = [
            domain if rank == "domain" else ranks.get(rank, "N/A")
            for rank in TAXONOMIC_RANKS
        ]
        return [sequence_id, aligned_clusters, aligned_clusters, *lineage, reason]

    @staticmethod
    def _write_raw(destination: Path, rows: list[list[str]]) -> None:
        header = [
            "sequence_name",
            "number_of_aligned_clusters",
            "number_of_genes_in_aligned_clusters",
            *TAXONOMIC_RANKS,
            "unassigned_reason",
        ]
        destination.write_text(
            "\n".join("\t".join(row) for row in [header, *rows]) + "\n",
            encoding="utf-8",
        )


if __name__ == "__main__":
    unittest.main()
