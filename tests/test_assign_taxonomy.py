import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "src" / "petra" / "workflow" / "scripts"
)
sys.path.insert(0, str(SCRIPTS_DIRECTORY))

from assign_taxonomy import TAXONOMIC_RANKS, assign_taxonomy  # noqa: E402


class AssignTaxonomyTest(unittest.TestCase):
    def test_assigns_from_half_of_aligned_clusters_and_reports_aligned_count(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            assembly = directory / "assembly.fasta"
            proteins = directory / "proteins.faa"
            cluster_sizes = directory / "cluster_sizes.tsv"
            alignment = directory / "lineages.tsv"
            destination = directory / "taxonomy.tsv"
            reasons = directory / "reasons.tsv"

            assembly.write_text(
                ">contig_one description\nACGT\n>contig_two\nACGT\n>no_genes\nACGT\n",
                encoding="utf-8",
            )
            proteins.write_text(
                ">contig_one_1\nM\n>contig_one_2\nM\n>contig_one_3\nM\n"
                ">contig_one_4\nM\n>contig_one_5\nM\n"
                ">contig_two_1\nM\n>contig_two_2\nM\n",
                encoding="utf-8",
            )
            self._write_cluster_sizes(
                cluster_sizes,
                {
                    "contig_one_1": 1,
                    "contig_one_2": 1,
                    "contig_one_3": 10,
                    "contig_one_4": 1,
                    "contig_one_5": 1,
                    "contig_two_1": 1,
                    "contig_two_2": 1,
                },
            )
            lineages = [
                self._lineage(domain="Bacteria", phylum="Pseudomonadota"),
                self._lineage(domain="Bacteria", phylum="Pseudomonadota"),
                self._lineage(domain="Archaea", phylum="Thermoproteota"),
                self._lineage(domain="N/A", phylum="N/A"),
                self._lineage(domain="Bacteria", phylum="Pseudomonadota"),
                self._lineage(domain="Archaea", phylum="Thermoproteota"),
            ]
            alignment.write_text(
                "".join(
                    f"{protein}\ttarget-{index}\t1e-20\t100\t{lineage}\n"
                    for index, (protein, lineage) in enumerate(
                        zip(
                            (
                                "contig_one_1",
                                "contig_one_2",
                                "contig_one_3",
                                "contig_one_4",
                                "contig_two_1",
                                "contig_two_2",
                            ),
                            lineages,
                            strict=True,
                        ),
                        start=1,
                    )
                ),
                encoding="utf-8",
            )

            assign_taxonomy(
                assembly, proteins, cluster_sizes, alignment, destination, reasons
            )

            rows = [line.split("\t") for line in destination.read_text().splitlines()]
            self.assertEqual(
                rows[0],
                [
                    "sequence_name",
                    "sequence_length",
                    "number_of_aligned_clusters",
                    "number_of_genes_in_aligned_clusters",
                    *TAXONOMIC_RANKS,
                ],
            )
            self.assertEqual(
                rows[1][0:6],
                ["contig_one", "4", "4", "13", "Bacteria", "N/A"],
            )
            self.assertEqual(rows[1][6], "Pseudomonadota")
            self.assertEqual(
                rows[2][0:6], ["contig_two", "4", "2", "2", "N/A", "N/A"]
            )
            self.assertEqual(rows[2][6], "N/A")
            self.assertEqual(
                rows[3],
                ["no_genes", "4", "0", "0", *(["N/A"] * len(TAXONOMIC_RANKS))],
            )
            self.assertEqual(
                [line.split("\t") for line in reasons.read_text().splitlines()],
                [
                    ["sequence_name", "unassigned_reason"],
                    ["contig_one", "-"],
                    ["contig_two", "no_consensus"],
                    ["no_genes", "-"],
                ],
            )

    def test_uses_one_assignment_per_protein_from_the_best_retained_hit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            assembly = directory / "assembly.fasta"
            proteins = directory / "proteins.faa"
            cluster_sizes = directory / "cluster_sizes.tsv"
            alignment = directory / "lineages.tsv"
            destination = directory / "taxonomy.tsv"
            reasons = directory / "reasons.tsv"

            assembly.write_text(">contig\nACGT\n", encoding="utf-8")
            proteins.write_text(">contig_1\nM\n>contig_2\nM\n", encoding="utf-8")
            self._write_cluster_sizes(cluster_sizes, {"contig_1": 1, "contig_2": 1})
            bacteria = self._lineage(domain="Bacteria")
            archaea = self._lineage(domain="Archaea")
            alignment.write_text(
                f"contig_1\tbest\t1e-30\t200\t{bacteria}\n"
                f"contig_1\tlower\t1e-20\t100\t{archaea}\n"
                f"contig_2\tbest\t1e-30\t200\t{bacteria}\n",
                encoding="utf-8",
            )

            assign_taxonomy(
                assembly, proteins, cluster_sizes, alignment, destination, reasons
            )

            row = destination.read_text(encoding="utf-8").splitlines()[1].split("\t")
            self.assertEqual(row[0:5], ["contig", "4", "2", "2", "Bacteria"])
            self.assertEqual(reasons.read_text().splitlines()[1], "contig\t-")

    def test_tied_best_hits_only_assign_ranks_on_which_they_agree(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            assembly = directory / "assembly.fasta"
            proteins = directory / "proteins.faa"
            cluster_sizes = directory / "cluster_sizes.tsv"
            alignment = directory / "lineages.tsv"
            destination = directory / "taxonomy.tsv"
            reasons = directory / "reasons.tsv"

            assembly.write_text(">contig\nACGT\n", encoding="utf-8")
            proteins.write_text(">contig_1\nM\n", encoding="utf-8")
            self._write_cluster_sizes(cluster_sizes, {"contig_1": 1})
            lineage_one = self._lineage(domain="Bacteria", phylum="Pseudomonadota")
            lineage_two = self._lineage(domain="Bacteria", phylum="Bacillota")
            alignment.write_text(
                f"contig_1\ttie-one\t1e-20\t100\t{lineage_one}\n"
                f"contig_1\ttie-two\t1e-20\t100\t{lineage_two}\n",
                encoding="utf-8",
            )

            assign_taxonomy(
                assembly, proteins, cluster_sizes, alignment, destination, reasons
            )

            row = destination.read_text(encoding="utf-8").splitlines()[1].split("\t")
            self.assertEqual(row[4], "Bacteria")
            self.assertEqual(row[6], "N/A")
            self.assertEqual(reasons.read_text().splitlines()[1], "contig\t-")

    def test_does_not_report_consensus_failure_without_contig_level_votes(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            assembly = directory / "assembly.fasta"
            proteins = directory / "proteins.faa"
            cluster_sizes = directory / "cluster_sizes.tsv"
            alignment = directory / "lineages.tsv"
            destination = directory / "taxonomy.tsv"
            reasons = directory / "reasons.tsv"

            assembly.write_text(">missing\nA\n>tied\nA\n", encoding="utf-8")
            proteins.write_text(">missing_1\nM\n>tied_1\nM\n", encoding="utf-8")
            self._write_cluster_sizes(cluster_sizes, {"missing_1": 1, "tied_1": 1})
            alignment.write_text(
                f"missing_1\tmissing\t1e-20\t100\t{self._lineage()}\n"
                f"tied_1\tbacteria\t1e-20\t100\t"
                f"{self._lineage(domain='Bacteria')}\n"
                f"tied_1\tarchaea\t1e-20\t100\t"
                f"{self._lineage(domain='Archaea')}\n",
                encoding="utf-8",
            )

            assign_taxonomy(
                assembly, proteins, cluster_sizes, alignment, destination, reasons
            )

            self.assertEqual(
                reasons.read_text().splitlines()[1:], ["missing\t-", "tied\t-"]
            )

    def test_requires_one_gene_in_aligned_clusters_per_configured_number_of_bases(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            assembly = directory / "assembly.fasta"
            proteins = directory / "proteins.faa"
            cluster_sizes = directory / "cluster_sizes.tsv"
            alignment = directory / "lineages.tsv"
            destination = directory / "taxonomy.tsv"
            reasons = directory / "reasons.tsv"

            assembly.write_text(
                ">at_threshold\nAAAA\nAAAA\n"
                ">rescued_by_members\nAAAAAAAAA\n"
                ">below_threshold\nAAAAAAAAA\n",
                encoding="utf-8",
            )
            proteins.write_text(
                ">at_threshold_1\nM\n>at_threshold_2\nM\n"
                ">rescued_by_members_1\nM\n>rescued_by_members_2\nM\n"
                ">below_threshold_1\nM\n>below_threshold_2\nM\n",
                encoding="utf-8",
            )
            self._write_cluster_sizes(
                cluster_sizes,
                {
                    "at_threshold_1": 1,
                    "at_threshold_2": 1,
                    "rescued_by_members_1": 2,
                    "rescued_by_members_2": 1,
                    "below_threshold_1": 1,
                    "below_threshold_2": 1,
                },
            )
            bacteria = self._lineage(domain="Bacteria")
            alignment.write_text(
                f"at_threshold_1\ttarget-1\t1e-20\t100\t{bacteria}\n"
                f"at_threshold_2\ttarget-2\t1e-20\t100\t{bacteria}\n"
                f"rescued_by_members_1\ttarget-3\t1e-20\t100\t{bacteria}\n"
                f"rescued_by_members_2\ttarget-4\t1e-20\t100\t{bacteria}\n"
                f"below_threshold_1\ttarget-5\t1e-20\t100\t{bacteria}\n"
                f"below_threshold_2\ttarget-6\t1e-20\t100\t{bacteria}\n",
                encoding="utf-8",
            )

            assign_taxonomy(
                assembly,
                proteins,
                cluster_sizes,
                alignment,
                destination,
                reasons,
                max_bases_per_aligned_gene=4,
            )

            rows = [line.split("\t") for line in destination.read_text().splitlines()]
            self.assertEqual(
                rows[1][0:5], ["at_threshold", "8", "2", "2", "Bacteria"]
            )
            self.assertEqual(
                rows[2][0:5], ["rescued_by_members", "9", "2", "3", "Bacteria"]
            )
            self.assertEqual(
                rows[3][0:5], ["below_threshold", "9", "2", "2", "N/A"]
            )
            self.assertEqual(
                reasons.read_text().splitlines()[1:],
                [
                    "at_threshold\t-",
                    "rescued_by_members\t-",
                    "below_threshold\tgene_density",
                ],
            )

    @staticmethod
    def _write_cluster_sizes(destination: Path, sizes: dict[str, int]) -> None:
        destination.write_text(
            "cluster_representative\tnumber_of_genes_in_cluster\n"
            + "".join(
                f"{representative_id}\t{size}\n"
                for representative_id, size in sizes.items()
            ),
            encoding="utf-8",
        )

    @staticmethod
    def _lineage(**ranks: str) -> str:
        return ";".join(ranks.get(rank, "N/A") for rank in TAXONOMIC_RANKS)


if __name__ == "__main__":
    unittest.main()
