import shutil
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "src" / "petra" / "workflow" / "scripts"
)
sys.path.insert(0, str(SCRIPTS_DIRECTORY))

from assign_taxonomy import TAXONOMIC_RANKS, assign_taxonomy  # noqa: E402

_DOMAIN_COLUMN = 3
_PHYLUM_COLUMN = _DOMAIN_COLUMN + TAXONOMIC_RANKS.index("phylum")


class AssignTaxonomyTest(unittest.TestCase):
    def test_assigns_ranks_supported_by_a_confident_majority(self) -> None:
        bacteria = self._lineage(domain="Bacteria", phylum="Pseudomonadota")
        archaea = self._lineage(domain="Archaea", phylum="Thermoproteota")
        rows, reasons = self._assign(
            {"contig": 4000, "no_genes": 100},
            self._protein_ids("contig", 10),
            self._hits("contig", [bacteria] * 8 + [archaea] * 2),
        )

        self.assertEqual(
            rows[0],
            [
                "sequence_name",
                "sequence_length",
                "number_of_aligned_clusters",
                *TAXONOMIC_RANKS,
            ],
        )
        self.assertEqual(rows[1][0:4], ["contig", "4000", "10", "Bacteria"])
        self.assertEqual(rows[1][_PHYLUM_COLUMN], "Pseudomonadota")
        self.assertEqual(
            rows[2], ["no_genes", "100", "0", *(["N/A"] * len(TAXONOMIC_RANKS))]
        )
        self.assertEqual(reasons[1:], [["contig", "-"], ["no_genes", "-"]])

    def test_unanimous_votes_need_a_minimum_number_of_clusters(self) -> None:
        bacteria = self._lineage(domain="Bacteria")
        cases = ((2, "N/A", "low_confidence"), (3, "Bacteria", "-"))
        for clusters, domain, reason in cases:
            with self.subTest(clusters=clusters):
                rows, reasons = self._assign(
                    {"contig": 100},
                    self._protein_ids("contig", clusters),
                    self._hits("contig", [bacteria] * clusters),
                )

                self.assertEqual(rows[1][_DOMAIN_COLUMN], domain)
                self.assertEqual(reasons[1], ["contig", reason])

    def test_reports_low_confidence_when_only_the_sample_size_rejects_the_majority(
        self,
    ) -> None:
        bacteria = self._lineage(domain="Bacteria")
        archaea = self._lineage(domain="Archaea")
        rows, reasons = self._assign(
            {"contig": 100},
            self._protein_ids("contig", 10),
            self._hits("contig", [bacteria] * 6 + [archaea] * 4),
        )

        self.assertEqual(rows[1][_DOMAIN_COLUMN], "N/A")
        self.assertEqual(reasons[1], ["contig", "low_confidence"])

    def test_reports_no_consensus_when_the_observed_share_rejects_the_majority(
        self,
    ) -> None:
        bacteria = self._lineage(domain="Bacteria")
        archaea = self._lineage(domain="Archaea")
        rows, reasons = self._assign(
            {"contig": 100},
            self._protein_ids("contig", 10),
            self._hits("contig", [bacteria] * 5 + [archaea] * 5),
        )

        self.assertEqual(rows[1][_DOMAIN_COLUMN], "N/A")
        self.assertEqual(reasons[1], ["contig", "no_consensus"])

    def test_alpha_sets_the_share_the_confidence_bound_must_reach(self) -> None:
        bacteria = self._lineage(domain="Bacteria")
        archaea = self._lineage(domain="Archaea")
        rows, reasons = self._assign(
            {"contig": 100},
            self._protein_ids("contig", 10),
            self._hits("contig", [bacteria] * 6 + [archaea] * 4),
            alpha=0.2,
        )

        self.assertEqual(rows[1][_DOMAIN_COLUMN], "Bacteria")
        self.assertEqual(reasons[1], ["contig", "-"])

    def test_rejects_an_alpha_outside_the_unit_interval(self) -> None:
        with self.assertRaisesRegex(ValueError, "alpha-majority threshold"):
            self._assign(
                {"contig": 100},
                self._protein_ids("contig", 1),
                self._hits("contig", [self._lineage(domain="Bacteria")]),
                alpha=0.0,
            )

    def test_uses_one_assignment_per_protein_from_the_best_retained_hit(self) -> None:
        bacteria = self._lineage(domain="Bacteria")
        archaea = self._lineage(domain="Archaea")
        hits = []
        for protein_id in self._protein_ids("contig", 3):
            hits.append((protein_id, "200", bacteria))
            hits.append((protein_id, "100", archaea))
        rows, reasons = self._assign(
            {"contig": 100}, self._protein_ids("contig", 3), hits
        )

        self.assertEqual(rows[1][0:4], ["contig", "100", "3", "Bacteria"])
        self.assertEqual(reasons[1], ["contig", "-"])

    def test_tied_best_hits_only_vote_on_ranks_on_which_they_agree(self) -> None:
        hits = []
        for protein_id in self._protein_ids("contig", 3):
            hits.append(
                (
                    protein_id,
                    "100",
                    self._lineage(domain="Bacteria", phylum="Bacillota"),
                )
            )
            hits.append(
                (
                    protein_id,
                    "100",
                    self._lineage(domain="Bacteria", phylum="Pseudomonadota"),
                )
            )
        rows, reasons = self._assign(
            {"contig": 100}, self._protein_ids("contig", 3), hits
        )

        self.assertEqual(rows[1][_DOMAIN_COLUMN], "Bacteria")
        self.assertEqual(rows[1][_PHYLUM_COLUMN], "N/A")
        self.assertEqual(reasons[1], ["contig", "-"])

    def test_does_not_report_a_failure_without_votes_for_any_rank(self) -> None:
        _, reasons = self._assign(
            {"missing": 100, "tied": 100},
            ["missing_1", "tied_1"],
            [
                ("missing_1", "100", self._lineage()),
                ("tied_1", "100", self._lineage(domain="Bacteria")),
                ("tied_1", "100", self._lineage(domain="Archaea")),
            ],
        )

        self.assertEqual(reasons[1:], [["missing", "-"], ["tied", "-"]])

    def _assign(
        self,
        sequence_lengths: dict[str, int],
        protein_ids: list[str],
        hits: list[tuple[str, str, str]],
        **options: float,
    ) -> tuple[list[list[str]], list[list[str]]]:
        directory = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, directory)
        assembly = directory / "assembly.fasta"
        proteins = directory / "proteins.faa"
        alignment = directory / "lineages.tsv"
        destination = directory / "taxonomy.tsv"
        reasons = directory / "reasons.tsv"

        assembly.write_text(
            "".join(
                f">{sequence_id}\n{'A' * length}\n"
                for sequence_id, length in sequence_lengths.items()
            ),
            encoding="utf-8",
        )
        proteins.write_text(
            "".join(f">{protein_id}\nM\n" for protein_id in protein_ids),
            encoding="utf-8",
        )
        alignment.write_text(
            "".join(
                f"{protein_id}\ttarget-{index}\t1e-20\t{bitscore}\t{lineage}\n"
                for index, (protein_id, bitscore, lineage) in enumerate(hits, start=1)
            ),
            encoding="utf-8",
        )

        assign_taxonomy(assembly, proteins, alignment, destination, reasons, **options)

        return (
            [line.split("\t") for line in destination.read_text().splitlines()],
            [line.split("\t") for line in reasons.read_text().splitlines()],
        )

    @staticmethod
    def _protein_ids(sequence_id: str, count: int) -> list[str]:
        return [f"{sequence_id}_{number}" for number in range(1, count + 1)]

    @staticmethod
    def _hits(sequence_id: str, lineages: list[str]) -> list[tuple[str, str, str]]:
        return [
            (f"{sequence_id}_{number}", "100", lineage)
            for number, lineage in enumerate(lineages, start=1)
        ]

    @staticmethod
    def _lineage(**ranks: str) -> str:
        return ";".join(ranks.get(rank, "N/A") for rank in TAXONOMIC_RANKS)


if __name__ == "__main__":
    unittest.main()
