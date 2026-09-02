import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "src" / "petra" / "workflow" / "scripts"
)
sys.path.insert(0, str(SCRIPTS_DIRECTORY))

from split_assembly_for_cmscan import split_assembly_for_cmscan  # noqa: E402


class SplitAssemblyForCmscanTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.assembly = self.directory / "assembly.fasta"

    def _chunks(self, count: int) -> list[Path]:
        return [self.directory / f"chunk.{index}.fasta" for index in range(count)]

    @staticmethod
    def _bases(chunk: Path) -> int:
        return sum(
            len(line.strip())
            for line in chunk.read_text(encoding="utf-8").splitlines()
            if not line.startswith(">")
        )

    def test_balances_chunks_by_bases_not_by_record_count(self) -> None:
        self.assembly.write_text(
            ">long\n" + "A" * 100 + "\n"
            ">medium-1\n" + "C" * 30 + "\n"
            ">medium-2\n" + "G" * 30 + "\n"
            ">medium-3\n" + "T" * 30 + "\n"
            ">medium-4\n" + "A" * 30 + "\n",
            encoding="utf-8",
        )
        chunks = self._chunks(2)

        split_assembly_for_cmscan(self.assembly, chunks)

        first, second = (chunk.read_text(encoding="utf-8") for chunk in chunks)
        self.assertIn(">long", first)
        self.assertEqual(first.count(">"), 1)
        self.assertEqual(second.count(">"), 4)
        self.assertEqual(
            [self._bases(chunk) for chunk in chunks],
            [100, 120],
        )

    def test_copies_records_verbatim_and_terminates_the_last_line(self) -> None:
        self.assembly.write_bytes(
            b">first description here\r\nAC\r\nGT\r\n>second\r\nTTTT"
        )
        chunks = self._chunks(2)

        split_assembly_for_cmscan(self.assembly, chunks)

        contents = sorted(chunk.read_bytes() for chunk in chunks)
        self.assertEqual(
            contents,
            [b">first description here\r\nAC\r\nGT\r\n", b">second\r\nTTTT\n"],
        )

    def test_creates_every_chunk_even_without_enough_sequences(self) -> None:
        self.assembly.write_text(">only\nACGT\n", encoding="utf-8")
        chunks = self._chunks(3)

        split_assembly_for_cmscan(self.assembly, chunks)

        self.assertTrue(all(chunk.is_file() for chunk in chunks))
        self.assertEqual(
            sorted(chunk.read_text(encoding="utf-8") for chunk in chunks),
            ["", "", ">only\nACGT\n"],
        )

    def test_assigns_the_same_way_on_every_run(self) -> None:
        self.assembly.write_text(
            "".join(f">sequence-{index}\n{'A' * (index + 1)}\n" for index in range(20)),
            encoding="utf-8",
        )
        first_chunks = self._chunks(4)
        split_assembly_for_cmscan(self.assembly, first_chunks)
        first = [chunk.read_text(encoding="utf-8") for chunk in first_chunks]

        second_chunks = [self.directory / f"repeat.{index}.fasta" for index in range(4)]
        split_assembly_for_cmscan(self.assembly, second_chunks)

        self.assertEqual(
            [chunk.read_text(encoding="utf-8") for chunk in second_chunks], first
        )

    def test_rejects_an_assembly_without_a_record(self) -> None:
        self.assembly.write_text("", encoding="utf-8")

        with self.assertRaisesRegex(ValueError, "No FASTA record"):
            split_assembly_for_cmscan(self.assembly, self._chunks(1))

    def test_rejects_sequence_data_before_the_first_header(self) -> None:
        self.assembly.write_text("ACGT\n>late\nACGT\n", encoding="utf-8")

        with self.assertRaisesRegex(ValueError, "before the first header"):
            split_assembly_for_cmscan(self.assembly, self._chunks(1))


if __name__ == "__main__":
    unittest.main()
