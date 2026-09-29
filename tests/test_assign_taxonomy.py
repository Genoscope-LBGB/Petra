import shutil
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "src" / "petra" / "workflow" / "scripts"
)
sys.path.insert(0, str(SCRIPTS_DIRECTORY))

from assign_taxonomy import (  # noqa: E402
    TAXONOMIC_RANKS,
    _minimum_aligned_unigenes,
    assign_taxonomy,
)

_DOMAIN_COLUMN = 4
_PHYLUM_COLUMN = _DOMAIN_COLUMN + TAXONOMIC_RANKS.index("phylum")


class AssignTaxonomyTest(unittest.TestCase):
    def test_assigns_ranks_supported_by_an_alpha_majority(self) -> None:
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
                "gc_percent",
                "number_of_aligned_unigenes",
                *TAXONOMIC_RANKS,
            ],
        )
        self.assertEqual(
            rows[1][0:5], ["contig", "4000", "0.00", "10", "Bacteria"]
        )
        self.assertEqual(rows[1][_PHYLUM_COLUMN], "Pseudomonadota")
        self.assertEqual(
            rows[2],
            ["no_genes", "100", "0.00", "0", *(["N/A"] * len(TAXONOMIC_RANKS))],
        )
        self.assertEqual(reasons[1:], [["contig", "-"], ["no_genes", "-"]])

    def test_reports_gc_percent_for_mixed_case_multiline_sequences(self) -> None:
        rows, _ = self._assign({"contig": "aCgT\nNNgc"}, [], [])

        self.assertEqual(rows[1][0:4], ["contig", "8", "50.00", "0"])

    def test_applies_the_minimum_and_length_anchor_to_the_whole_sequence(self) -> None:
        bacteria = self._lineage(domain="Bacteria", phylum="Pseudomonadota")
        cases = (
            (100, 2, "N/A", "insufficient_support"),
            (50_000, 3, "Bacteria", "-"),
            (50_001, 3, "N/A", "insufficient_support"),
            (50_001, 4, "Bacteria", "-"),
        )
        for length, unigenes, domain, reason in cases:
            with self.subTest(length=length, unigenes=unigenes):
                rows, reasons = self._assign(
                    {"contig": length},
                    self._protein_ids("contig", unigenes),
                    self._hits("contig", [bacteria] * unigenes),
                )
                self.assertEqual(rows[1][_DOMAIN_COLUMN], domain)
                self.assertEqual(reasons[1], ["contig", reason])
                if domain == "N/A":
                    self.assertEqual(
                        rows[1][_DOMAIN_COLUMN:], ["N/A"] * len(TAXONOMIC_RANKS)
                    )

    def test_logarithmic_support_interpolates_and_rounds_up(self) -> None:
        cases = (
            (0, 3),
            (100_000, 4),
            (500_000, 6),
            (500_001, 7),
            (5_000_000, 9),
            (50_000_000, 12),
            (50_000_001, 13),
        )
        for length, minimum in cases:
            with self.subTest(length=length):
                self.assertEqual(
                    _minimum_aligned_unigenes(length, 3, 50_000, 3), minimum
                )

    def test_support_counts_aligned_unigenes_not_predictions_or_hits(self) -> None:
        bacteria = self._lineage(domain="Bacteria")
        hits = self._hits("contig", [bacteria] * 3)
        rows, reasons = self._assign(
            {"contig": 50_001},
            self._protein_ids("contig", 10),
            hits + hits,
        )
        self.assertEqual(rows[1][3:5], ["3", "N/A"])
        self.assertEqual(reasons[1], ["contig", "insufficient_support"])

    def test_support_is_not_required_again_for_each_rank(self) -> None:
        rows, reasons = self._assign(
            {"contig": 500_000},
            self._protein_ids("contig", 7),
            self._hits(
                "contig",
                [self._lineage(domain="Bacteria", phylum="Pseudomonadota")] * 4
                + [self._lineage(domain="Bacteria", phylum="Bacillota")]
                + [self._lineage(domain="Archaea", phylum="Thermoproteota")] * 2,
            ),
        )
        self.assertEqual(rows[1][_DOMAIN_COLUMN], "Bacteria")
        self.assertEqual(rows[1][_PHYLUM_COLUMN], "Pseudomonadota")
        self.assertEqual(reasons[1], ["contig", "-"])

    def test_accepts_a_narrow_observed_majority(self) -> None:
        rows, reasons = self._assign(
            {"contig": 100},
            self._protein_ids("contig", 100),
            self._hits(
                "contig",
                [self._lineage(domain="Bacteria")] * 51
                + [self._lineage(domain="Archaea")] * 49,
            ),
        )
        self.assertEqual(rows[1][_DOMAIN_COLUMN], "Bacteria")
        self.assertEqual(reasons[1], ["contig", "-"])

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

    def test_leaves_ranks_below_an_unresolved_rank_unassigned(self) -> None:
        virus = self._lineage(domain="Viruses", **{"class": "Naldaviricetes"})
        animal = self._lineage(
            domain="Eukaryota", kingdom="Metazoa", phylum="Arthropoda"
        )
        rows, reasons = self._assign(
            {"contig": 100},
            self._protein_ids("contig", 6),
            self._hits("contig", [virus] * 3 + [animal] * 3),
        )

        self.assertEqual(rows[1][_DOMAIN_COLUMN:], ["N/A"] * len(TAXONOMIC_RANKS))
        self.assertEqual(reasons[1], ["contig", "no_consensus"])

    def test_descends_past_a_rank_absent_from_the_winning_clade(self) -> None:
        ciliate = self._lineage(domain="Eukaryota", phylum="Ciliophora")
        animal = self._lineage(
            domain="Eukaryota", kingdom="Metazoa", phylum="Cnidaria"
        )
        rows, _ = self._assign(
            {"contig": 100},
            self._protein_ids("contig", 4),
            self._hits("contig", [ciliate] * 3 + [animal]),
        )

        self.assertEqual(
            rows[1][_DOMAIN_COLUMN : _PHYLUM_COLUMN + 1],
            ["Eukaryota", "N/A", "Ciliophora"],
        )

    def test_a_share_equal_to_alpha_is_not_a_majority(self) -> None:
        bacteria = self._lineage(domain="Bacteria")
        for options, domain in (({}, "N/A"), ({"alpha": 0.49}, "Bacteria")):
            with self.subTest(options=options):
                rows, reasons = self._assign(
                    {"contig": 100},
                    self._protein_ids("contig", 10),
                    self._hits("contig", [bacteria] * 5 + [self._lineage()] * 5),
                    **options,
                )
                self.assertEqual(rows[1][_DOMAIN_COLUMN], domain)
                self.assertEqual(
                    reasons[1],
                    ["contig", "-" if domain != "N/A" else "no_consensus"],
                )

    def test_custom_support_curve_controls_assignment(self) -> None:
        # At 1 kb require four unigenes; at 10 kb require six.
        bacteria = self._lineage(domain="Bacteria")
        for length, unigenes, domain in (
            (1_000, 3, "N/A"),
            (1_000, 4, "Bacteria"),
            (10_000, 5, "N/A"),
            (10_000, 6, "Bacteria"),
        ):
            with self.subTest(length=length, unigenes=unigenes):
                rows, _ = self._assign(
                    {"contig": length},
                    self._protein_ids("contig", unigenes),
                    self._hits("contig", [bacteria] * unigenes),
                    min_aligned_unigenes=4,
                    unigene_support_length=1_000,
                    unigenes_per_decade=2,
                )
                self.assertEqual(rows[1][_DOMAIN_COLUMN], domain)

    def test_rejects_invalid_support_parameters(self) -> None:
        for options in (
            {"min_aligned_unigenes": 2},
            {"unigene_support_length": 0},
            {"unigenes_per_decade": 0},
        ):
            with self.subTest(options=options), self.assertRaises(ValueError):
                self._assign({"contig": 100}, [], [], **options)

    def test_rejects_an_alpha_outside_the_open_unit_interval(self) -> None:
        for alpha in (0.0, 1.0):
            with (
                self.subTest(alpha=alpha),
                self.assertRaisesRegex(ValueError, "alpha-majority threshold"),
            ):
                self._assign(
                    {"contig": 100},
                    self._protein_ids("contig", 1),
                    self._hits("contig", [self._lineage(domain="Bacteria")]),
                    alpha=alpha,
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

        self.assertEqual(
            rows[1][0:5], ["contig", "100", "0.00", "3", "Bacteria"]
        )
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

    def test_reports_each_aligned_gene_under_its_assembly_sequence(self) -> None:
        bacteria = self._lineage(domain="Bacteria", phylum="Pseudomonadota")
        archaea = self._lineage(domain="Archaea")
        self._assign(
            {"contig_1": 100, "empty": 100},
            ["contig_1_1", "contig_1_2", "contig_1_3"],
            [
                (
                    "contig_1_2",
                    "100",
                    self._lineage(domain="Bacteria", phylum="Bacillota"),
                ),
                ("contig_1_1", "50", archaea),
                ("contig_1_2", "100", bacteria),
                ("contig_1_1", "100", bacteria),
            ],
        )

        self.assertEqual(
            self.gene_taxonomy_rows,
            [
                ["contig_1"],
                ["", "contig_1_1", "Bacteria;Pseudomonadota"],
                ["", "contig_1_2", "Bacteria"],
                ["empty"],
            ],
        )

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
        self.assertEqual(
            self.gene_taxonomy_rows,
            [
                ["missing"],
                ["", "missing_1", "N/A"],
                ["tied"],
                ["", "tied_1", "N/A"],
            ],
        )

    def _assign(
        self,
        sequences: dict[str, int | str],
        protein_ids: list[str],
        hits: list[tuple[str, str, str]],
        **options: int | float,
    ) -> tuple[list[list[str]], list[list[str]]]:
        directory = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, directory)
        assembly = directory / "assembly.fasta"
        proteins = directory / "proteins.faa"
        alignment = directory / "lineages.tsv"
        destination = directory / "taxonomy.tsv"
        reasons = directory / "reasons.tsv"
        gene_taxonomy = directory / "gene-taxonomy.tsv"

        assembly.write_text(
            "".join(
                f">{sequence_id}\n"
                f"{'A' * sequence if isinstance(sequence, int) else sequence}\n"
                for sequence_id, sequence in sequences.items()
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

        assign_taxonomy(
            assembly,
            proteins,
            alignment,
            destination,
            reasons,
            gene_taxonomy,
            **options,
        )
        self.gene_taxonomy_rows = [
            line.split("\t") for line in gene_taxonomy.read_text().splitlines()
        ]

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
