import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "src" / "petra" / "workflow" / "scripts"
)
sys.path.insert(0, str(SCRIPTS_DIRECTORY))

from extract_rrna_genes import extract_rrna_genes, parse_tblout  # noqa: E402

_HEADER = "#idx target name accession query name accession clan name\n"


def _row(
    *,
    family: str,
    accession: str,
    sequence: str,
    start: int,
    end: int,
    strand: str,
    score: str = "1500.0",
    evalue: str = "1e-20",
    overlap: str = "*",
    clan: str = "CL00111",
) -> str:
    """Build one cmscan '--fmt 2' tblout row, 27 whitespace-separated fields."""
    fields = [
        "1",
        family,
        accession,
        sequence,
        "-",
        clan,
        "cm",
        "1",
        "1800",
        str(start),
        str(end),
        strand,
        "no",
        "1",
        "0.50",
        "0.0",
        score,
        evalue,
        "!",
        overlap,
        "-",
        "-",
        "-",
        "-",
        "-",
        "-",
        "description",
    ]
    return " ".join(fields) + "\n"


class ParseTbloutTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.tblout = self.directory / "predictions.rrna.tblout"

    def test_drops_clan_competitors_and_orients_coordinates(self) -> None:
        self.tblout.write_text(
            _HEADER
            + _row(
                family="SSU_rRNA_eukarya",
                accession="RF01960",
                sequence="contig-1",
                start=200,
                end=51,
                strand="-",
            )
            + _row(
                family="SSU_rRNA_bacteria",
                accession="RF00177",
                sequence="contig-1",
                start=60,
                end=190,
                strand="+",
                score="900.0",
                overlap="=",
            ),
            encoding="utf-8",
        )

        genes = parse_tblout(self.tblout)

        self.assertEqual(len(genes), 1)
        gene = genes[0]
        self.assertEqual((gene.start, gene.end, gene.strand), (51, 200, "-"))
        self.assertEqual(gene.subunit, "SSU")
        self.assertEqual(gene.gene_id, "rrna_1")

    def test_classifies_the_subunit_of_every_covariance_model(self) -> None:
        rows = [
            ("RF00001", "5S_rRNA", None),
            ("RF00002", "5_8S_rRNA", None),
            ("RF00177", "SSU_rRNA_bacteria", "SSU"),
            ("RF01959", "SSU_rRNA_archaea", "SSU"),
            ("RF01960", "SSU_rRNA_eukarya", "SSU"),
            ("RF02542", "SSU_rRNA_microsporidia", "SSU"),
            ("RF02545", "SSU_trypano_mito", "SSU"),
            ("RF02540", "LSU_rRNA_archaea", "LSU"),
            ("RF02541", "LSU_rRNA_bacteria", "LSU"),
            ("RF02543", "LSU_rRNA_eukarya", "LSU"),
            ("RF02546", "LSU_trypano_mito", "LSU"),
        ]
        self.tblout.write_text(
            "".join(
                _row(
                    family=family,
                    accession=accession,
                    sequence="contig-1",
                    start=index * 10 + 1,
                    end=index * 10 + 9,
                    strand="+",
                )
                for index, (accession, family, _) in enumerate(rows)
            ),
            encoding="utf-8",
        )

        genes = parse_tblout(self.tblout)

        self.assertEqual(
            [gene.subunit for gene in genes], [subunit for _, _, subunit in rows]
        )

    def test_rejects_a_table_written_without_fmt_2(self) -> None:
        self.tblout.write_text(
            "SSU_rRNA_eukarya RF01960 contig-1 - cm 1 1800 51 200 + no 1 0.50 "
            "0.0 1500.0 1e-20 ! description\n",
            encoding="utf-8",
        )

        with self.assertRaisesRegex(ValueError, r"--fmt 2"):
            parse_tblout(self.tblout)


class ExtractRnaGenesTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.tblout = self.directory / "predictions.rrna.tblout"
        self.assembly = self.directory / "assembly.fasta"
        self.outputs = {
            "gff": self.directory / "predictions.rrna.gff",
            "genes_table": self.directory / "predictions.rrna.genes.tsv",
            "ssu": self.directory / "predictions.rrna.ssu.fasta",
            "lsu": self.directory / "predictions.rrna.lsu.fasta",
        }

    def _extract(self) -> None:
        extract_rrna_genes(self.assembly, self.tblout, **self.outputs)

    def test_extracts_both_strands_at_the_reported_coordinates(self) -> None:
        self.assembly.write_text(
            ">contig-1 first\nAAACCCGGGTTT\n>contig-2\nACGTNNNNACGT\n",
            encoding="utf-8",
        )
        self.tblout.write_text(
            _row(
                family="SSU_rRNA_eukarya",
                accession="RF01960",
                sequence="contig-1",
                start=4,
                end=9,
                strand="+",
            )
            + _row(
                family="LSU_rRNA_eukarya",
                accession="RF02543",
                sequence="contig-2",
                start=8,
                end=1,
                strand="-",
            ),
            encoding="utf-8",
        )

        self._extract()

        self.assertEqual(
            self.outputs["ssu"].read_text(encoding="utf-8"), ">rrna_1\nCCCGGG\n"
        )
        self.assertEqual(
            self.outputs["lsu"].read_text(encoding="utf-8"), ">rrna_2\nNNNNACGT\n"
        )

    def test_reports_unclassifiable_families_in_the_gff_only(self) -> None:
        self.assembly.write_text(">contig-1\n" + "A" * 200 + "\n", encoding="utf-8")
        self.tblout.write_text(
            _row(
                family="5S_rRNA",
                accession="RF00001",
                sequence="contig-1",
                start=10,
                end=120,
                strand="+",
            ),
            encoding="utf-8",
        )

        self._extract()

        gff = self.outputs["gff"].read_text(encoding="utf-8").splitlines()
        self.assertEqual(gff[0], "##gff-version 3")
        self.assertEqual(
            gff[1].split("\t")[:8],
            ["contig-1", "cmscan", "rRNA", "10", "120", "1500.0", "+", "."],
        )
        self.assertIn("rfam_accession=RF00001", gff[1])
        self.assertEqual(
            self.outputs["genes_table"].read_text(encoding="utf-8").splitlines(),
            [
                (
                    "gene_id\tsequence_name\tgene_start\tgene_end\tstrand\t"
                    "rfam_accession\tsubunit\tcmscan_score"
                )
            ],
        )
        self.assertEqual(self.outputs["ssu"].read_text(encoding="utf-8"), "")
        self.assertEqual(self.outputs["lsu"].read_text(encoding="utf-8"), "")

    def test_writes_every_output_when_no_gene_was_detected(self) -> None:
        self.assembly.write_text(">contig-1\nACGT\n", encoding="utf-8")
        self.tblout.write_text(_HEADER, encoding="utf-8")

        self._extract()

        self.assertEqual(
            self.outputs["gff"].read_text(encoding="utf-8"), "##gff-version 3\n"
        )
        self.assertEqual(self.outputs["ssu"].read_text(encoding="utf-8"), "")

    def test_rejects_a_hit_beyond_the_end_of_its_sequence(self) -> None:
        self.assembly.write_text(">contig-1\nACGT\n", encoding="utf-8")
        self.tblout.write_text(
            _row(
                family="SSU_rRNA_eukarya",
                accession="RF01960",
                sequence="contig-1",
                start=1,
                end=10,
                strand="+",
            ),
            encoding="utf-8",
        )

        with self.assertRaisesRegex(ValueError, "ends at 10"):
            self._extract()

    def test_rejects_a_hit_on_a_sequence_absent_from_the_assembly(self) -> None:
        self.assembly.write_text(">contig-1\nACGT\n", encoding="utf-8")
        self.tblout.write_text(
            _row(
                family="SSU_rRNA_eukarya",
                accession="RF01960",
                sequence="contig-9",
                start=1,
                end=4,
                strand="+",
            ),
            encoding="utf-8",
        )

        with self.assertRaisesRegex(ValueError, "absent from"):
            self._extract()


if __name__ == "__main__":
    unittest.main()
