import csv
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "src" / "petra" / "workflow" / "scripts"
)
sys.path.insert(0, str(SCRIPTS_DIRECTORY))

from assign_taxonomy import TAXONOMIC_RANKS  # noqa: E402
from split_assembly_by_domain import split_assembly_by_domain  # noqa: E402


class SplitAssemblyByDomainTest(unittest.TestCase):
    def test_partitions_records_verbatim_in_assembly_order(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            assembly = directory / "assembly.fasta"
            taxonomy = directory / "taxonomy.tsv"
            outputs = self._outputs(directory)
            assembly.write_bytes(
                b">virus first\r\nAC\r\n"
                b">eukaryote description\r\naa\r\n"
                b">viroid\r\nGT\r\n"
                b">bacterium\r\nCC\r\n"
                b">archaeon\r\nTA\r\n"
                b">unassigned\r\nNN"
            )
            self._write_taxonomy(
                taxonomy,
                [
                    ("unassigned", "N/A", "whokaryote"),
                    ("archaeon", "Archaea", "-"),
                    ("bacterium", "Bacteria", "-"),
                    ("viroid", "Viroids", "-"),
                    ("eukaryote", "Eukaryota", "-"),
                    ("virus", "Viruses", "-"),
                ],
            )

            split_assembly_by_domain(assembly, taxonomy, **outputs)

            self.assertEqual(
                outputs["eukaryota"].read_bytes(),
                b">eukaryote description\r\naa\r\n",
            )
            self.assertEqual(outputs["bacteria"].read_bytes(), b">bacterium\r\nCC\r\n")
            self.assertEqual(
                outputs["viruses"].read_bytes(),
                b">virus first\r\nAC\r\n>viroid\r\nGT\r\n",
            )
            self.assertEqual(outputs["archaea"].read_bytes(), b">archaeon\r\nTA\r\n")
            self.assertEqual(outputs["unassigned"].read_bytes(), b">unassigned\r\nNN")

    def test_creates_empty_files_for_unrepresented_domains(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            assembly = directory / "assembly.fasta"
            taxonomy = directory / "taxonomy.tsv"
            outputs = self._outputs(directory)
            assembly.write_text(">bacterium\nAC\n", encoding="utf-8")
            self._write_taxonomy(taxonomy, [("bacterium", "Bacteria", "-")])

            split_assembly_by_domain(assembly, taxonomy, **outputs)

            self.assertEqual(outputs["bacteria"].read_text(), ">bacterium\nAC\n")
            for output_name in ("eukaryota", "viruses", "archaea", "unassigned"):
                self.assertEqual(outputs[output_name].read_bytes(), b"")

    def test_rejects_assembly_sequence_absent_from_taxonomy(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            assembly = directory / "assembly.fasta"
            taxonomy = directory / "taxonomy.tsv"
            assembly.write_text(">missing\nAC\n", encoding="utf-8")
            self._write_taxonomy(taxonomy, [])

            with self.assertRaisesRegex(ValueError, "absent from the final taxonomy"):
                split_assembly_by_domain(assembly, taxonomy, **self._outputs(directory))

    def test_rejects_taxonomy_sequence_absent_from_assembly(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            assembly = directory / "assembly.fasta"
            taxonomy = directory / "taxonomy.tsv"
            assembly.write_text(">present\nAC\n", encoding="utf-8")
            self._write_taxonomy(
                taxonomy,
                [("present", "Bacteria", "-"), ("missing", "Archaea", "-")],
            )

            with self.assertRaisesRegex(ValueError, "absent from the assembly"):
                split_assembly_by_domain(assembly, taxonomy, **self._outputs(directory))

    @staticmethod
    def _outputs(directory: Path) -> dict[str, Path]:
        return {
            name: directory / f"{name}.fasta"
            for name in (
                "eukaryota",
                "bacteria",
                "viruses",
                "archaea",
                "unassigned",
            )
        }

    @staticmethod
    def _write_taxonomy(destination: Path, rows: list[tuple[str, str, str]]) -> None:
        header = [
            "sequence_name",
            "number_of_aligned_clusters",
            "number_of_genes_in_aligned_clusters",
            *TAXONOMIC_RANKS,
            "unassigned_reason",
        ]
        with destination.open("w", encoding="utf-8", newline="") as output_handle:
            writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            for sequence_id, domain, reason in rows:
                writer.writerow(
                    [
                        sequence_id,
                        "1",
                        "1",
                        domain,
                        *(["N/A"] * (len(TAXONOMIC_RANKS) - 1)),
                        reason,
                    ]
                )


if __name__ == "__main__":
    unittest.main()
