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
            alignment = directory / "lineages.tsv"
            destination = directory / "taxonomy.tsv"

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

            assign_taxonomy(assembly, proteins, alignment, destination)

            rows = [line.split("\t") for line in destination.read_text().splitlines()]
            self.assertEqual(
                rows[0],
                ["sequence_name", "number_of_aligned_clusters", *TAXONOMIC_RANKS],
            )
            self.assertEqual(rows[1][0:4], ["contig_one", "4", "Bacteria", "N/A"])
            self.assertEqual(rows[1][4], "Pseudomonadota")
            self.assertEqual(rows[2][0:4], ["contig_two", "2", "N/A", "N/A"])
            self.assertEqual(rows[2][4], "N/A")
            self.assertEqual(rows[3], ["no_genes", "0", *(["N/A"] * 15)])

    def test_uses_one_assignment_per_protein_from_the_best_retained_hit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            assembly = directory / "assembly.fasta"
            proteins = directory / "proteins.faa"
            alignment = directory / "lineages.tsv"
            destination = directory / "taxonomy.tsv"

            assembly.write_text(">contig\nACGT\n", encoding="utf-8")
            proteins.write_text(">contig_1\nM\n>contig_2\nM\n", encoding="utf-8")
            bacteria = self._lineage(domain="Bacteria")
            archaea = self._lineage(domain="Archaea")
            alignment.write_text(
                f"contig_1\tbest\t1e-30\t200\t{bacteria}\n"
                f"contig_1\tlower\t1e-20\t100\t{archaea}\n"
                f"contig_2\tbest\t1e-30\t200\t{bacteria}\n",
                encoding="utf-8",
            )

            assign_taxonomy(assembly, proteins, alignment, destination)

            row = destination.read_text(encoding="utf-8").splitlines()[1].split("\t")
            self.assertEqual(row[0:3], ["contig", "2", "Bacteria"])

    def test_tied_best_hits_only_assign_ranks_on_which_they_agree(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            assembly = directory / "assembly.fasta"
            proteins = directory / "proteins.faa"
            alignment = directory / "lineages.tsv"
            destination = directory / "taxonomy.tsv"

            assembly.write_text(">contig\nACGT\n", encoding="utf-8")
            proteins.write_text(">contig_1\nM\n", encoding="utf-8")
            lineage_one = self._lineage(domain="Bacteria", phylum="Pseudomonadota")
            lineage_two = self._lineage(domain="Bacteria", phylum="Bacillota")
            alignment.write_text(
                f"contig_1\ttie-one\t1e-20\t100\t{lineage_one}\n"
                f"contig_1\ttie-two\t1e-20\t100\t{lineage_two}\n",
                encoding="utf-8",
            )

            assign_taxonomy(assembly, proteins, alignment, destination)

            row = destination.read_text(encoding="utf-8").splitlines()[1].split("\t")
            self.assertEqual(row[2], "Bacteria")
            self.assertEqual(row[4], "N/A")

    def test_requires_one_aligned_cluster_per_configured_number_of_bases(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            assembly = directory / "assembly.fasta"
            proteins = directory / "proteins.faa"
            alignment = directory / "lineages.tsv"
            destination = directory / "taxonomy.tsv"

            assembly.write_text(
                ">at_threshold\nAAAA\nAAAA\n>below_threshold\nAAAAAAAAA\n",
                encoding="utf-8",
            )
            proteins.write_text(
                ">at_threshold_1\nM\n>at_threshold_2\nM\n"
                ">below_threshold_1\nM\n>below_threshold_2\nM\n",
                encoding="utf-8",
            )
            bacteria = self._lineage(domain="Bacteria")
            alignment.write_text(
                f"at_threshold_1\ttarget-1\t1e-20\t100\t{bacteria}\n"
                f"at_threshold_2\ttarget-2\t1e-20\t100\t{bacteria}\n"
                f"below_threshold_1\ttarget-3\t1e-20\t100\t{bacteria}\n"
                f"below_threshold_2\ttarget-4\t1e-20\t100\t{bacteria}\n",
                encoding="utf-8",
            )

            assign_taxonomy(
                assembly,
                proteins,
                alignment,
                destination,
                max_bases_per_aligned_gene=4,
            )

            rows = [line.split("\t") for line in destination.read_text().splitlines()]
            self.assertEqual(rows[1][0:3], ["at_threshold", "2", "Bacteria"])
            self.assertEqual(rows[2][0:3], ["below_threshold", "2", "N/A"])

    @staticmethod
    def _lineage(**ranks: str) -> str:
        return ";".join(ranks.get(rank, "N/A") for rank in TAXONOMIC_RANKS)


if __name__ == "__main__":
    unittest.main()
