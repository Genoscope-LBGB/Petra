import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "src" / "petra" / "workflow" / "scripts"
)
sys.path.insert(0, str(SCRIPTS_DIRECTORY))

from select_rrna_best_hits import select_rrna_best_hits  # noqa: E402

_GENE_HEADER = (
    "gene_id\tsequence_name\tgene_start\tgene_end\tstrand\trfam_accession\t"
    "subunit\tcmscan_score\n"
)


def _hit(
    *,
    query: str,
    subject: str,
    identity: str = "99.000",
    alignment_length: int = 1000,
    bitscore: str = "1800",
    query_coverage: str = "98",
    taxid: str = "9606",
) -> str:
    """Build one 'blastn -outfmt "6 std qlen slen qcovs staxid"' row."""
    fields = [
        query,
        subject,
        identity,
        str(alignment_length),
        "5",
        "1",
        "1",
        str(alignment_length),
        "1",
        str(alignment_length),
        "0.0",
        bitscore,
        "1800",
        "1800",
        query_coverage,
        taxid,
    ]
    return "\t".join(fields) + "\n"


class SelectRnaBestHitsTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.genes = self.directory / "predictions.rrna.genes.tsv"
        self.ssu = self.directory / "ssu.blast.tsv"
        self.lsu = self.directory / "lsu.blast.tsv"
        self.destination = self.directory / "predictions.rrna.hits.tsv"
        self.ssu.write_text("", encoding="utf-8")
        self.lsu.write_text("", encoding="utf-8")

    def _write_genes(self, *rows: str) -> None:
        self.genes.write_text(_GENE_HEADER + "".join(rows), encoding="utf-8")

    def _select(self, minimum: int = 200) -> list[list[str]]:
        select_rrna_best_hits(
            self.genes,
            self.ssu,
            self.lsu,
            self.destination,
            minimum_alignment_length=minimum,
        )
        return [
            line.split("\t")
            for line in self.destination.read_text(encoding="utf-8").splitlines()
        ]

    def test_breaks_bitscore_ties_on_identity_length_then_accession(self) -> None:
        self._write_genes("rrna_1\tcontig-1\t1\t1800\t+\tRF01960\tSSU\t1500.0\n")
        self.ssu.write_text(
            _hit(query="rrna_1", subject="zzz-first-in-file", identity="98.000")
            + _hit(query="rrna_1", subject="bbb-same-identity", identity="99.000")
            + _hit(
                query="rrna_1",
                subject="aaa-same-identity",
                identity="99.000",
                alignment_length=900,
            ),
            encoding="utf-8",
        )

        rows = self._select()

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0][8], "bbb-same-identity")

    def test_prefers_the_highest_bitscore_over_identity(self) -> None:
        self._write_genes("rrna_1\tcontig-1\t1\t1800\t+\tRF01960\tSSU\t1500.0\n")
        self.ssu.write_text(
            _hit(
                query="rrna_1",
                subject="short-perfect",
                identity="100.000",
                alignment_length=400,
                bitscore="700",
            )
            + _hit(
                query="rrna_1", subject="long-good", identity="97.000", bitscore="1750"
            ),
            encoding="utf-8",
        )

        rows = self._select()

        self.assertEqual(rows[0][8], "long-good")

    def test_discards_alignments_below_the_minimum_length(self) -> None:
        self._write_genes(
            "rrna_1\tcontig-1\t1\t1800\t+\tRF01960\tSSU\t1500.0\n"
            "rrna_2\tcontig-2\t1\t3500\t+\tRF02543\tLSU\t2500.0\n"
        )
        self.ssu.write_text(
            _hit(query="rrna_1", subject="too-short", alignment_length=199),
            encoding="utf-8",
        )
        self.lsu.write_text(
            _hit(query="rrna_2", subject="at-the-limit", alignment_length=200),
            encoding="utf-8",
        )

        rows = self._select()

        self.assertEqual([row[0] for row in rows], ["rrna_2"])
        self.assertEqual(rows[0][8], "at-the-limit")

    def test_keeps_each_subunit_inside_its_own_database(self) -> None:
        self._write_genes("rrna_1\tcontig-1\t1\t1800\t+\tRF01960\tSSU\t1500.0\n")
        self.lsu.write_text(
            _hit(query="rrna_1", subject="wrong-database"), encoding="utf-8"
        )

        self.assertEqual(self._select(), [])

    def test_normalizes_a_missing_or_multiple_taxid(self) -> None:
        self._write_genes(
            "rrna_1\tcontig-1\t1\t1800\t+\tRF01960\tSSU\t1500.0\n"
            "rrna_2\tcontig-2\t1\t3500\t+\tRF02543\tLSU\t2500.0\n"
        )
        self.ssu.write_text(
            _hit(query="rrna_1", subject="without-taxid", taxid=""), encoding="utf-8"
        )
        self.lsu.write_text(
            _hit(query="rrna_2", subject="several-taxids", taxid="9606;9605"),
            encoding="utf-8",
        )

        rows = self._select()

        self.assertEqual([row[-1] for row in rows], ["0", "9606"])

    def test_writes_nothing_when_no_gene_was_classified(self) -> None:
        self._write_genes()

        self.assertEqual(self._select(), [])

    def test_rejects_an_alignment_written_without_the_extra_columns(self) -> None:
        self._write_genes("rrna_1\tcontig-1\t1\t1800\t+\tRF01960\tSSU\t1500.0\n")
        self.ssu.write_text(
            "rrna_1\tsubject\t99.0\t1000\t5\t1\t1\t1000\t1\t1000\t0.0\t1800\n",
            encoding="utf-8",
        )

        with self.assertRaisesRegex(ValueError, "qcovs staxid"):
            self._select()

    def test_rejects_an_unexpected_gene_table(self) -> None:
        self.genes.write_text("gene_id\tsubunit\nrrna_1\tSSU\n", encoding="utf-8")

        with self.assertRaisesRegex(ValueError, "Unexpected gene table columns"):
            self._select()


if __name__ == "__main__":
    unittest.main()
