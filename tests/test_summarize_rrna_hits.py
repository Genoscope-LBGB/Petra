import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "src" / "petra" / "workflow" / "scripts"
)
sys.path.insert(0, str(SCRIPTS_DIRECTORY))

from summarize_rrna_hits import summarize_rrna_hits

_HEADER = (
    "sequence_name\tgene_start\tgene_end\tstrand\trfam_accession\tsubunit\t"
    "cmscan_score\tsubject\tpercent_identity\talignment_length\t"
    "query_coverage\tlineage\tnon_host"
)


def _row(
    *,
    gene_id: str = "rrna_1",
    sequence: str = "contig-1",
    start: int = 100,
    end: int = 1800,
    subunit: str = "SSU",
    subject: str = "AB000001",
    taxid: str = "8022",
    lineage: str = (
        "Eukaryota;Metazoa;Chordata;N/A;N/A;Actinopteri;N/A;Salmoniformes;N/A;"
        "Salmonidae;N/A;Oncorhynchus;N/A;Oncorhynchus mykiss;N/A"
    ),
) -> str:
    fields = [
        gene_id,
        sequence,
        str(start),
        str(end),
        "+",
        "RF01960",
        subunit,
        "1500.0",
        subject,
        "99.520",
        "1750",
        "98",
        taxid,
        lineage,
    ]
    return "\t".join(fields) + "\n"


class SummarizeRnaHitsTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.lineages = self.directory / "predictions.rrna.lineages.tsv"
        self.destination = self.directory / "predictions.rrna.tsv"
        self.warnings = self.directory / "predictions.rrna.warnings.txt"

    def _summarize(
        self, host_phylum: str | None, host_taxon: str | None = None
    ) -> list[str]:
        summarize_rrna_hits(
            self.lineages,
            self.destination,
            self.warnings,
            host_taxon=host_taxon,
            host_phylum=host_phylum,
        )
        return self.destination.read_text(encoding="utf-8").splitlines()

    def test_reports_a_pipe_joined_lineage_without_a_host(self) -> None:
        self.lineages.write_text(_row(), encoding="utf-8")

        lines = self._summarize(None)

        self.assertEqual(lines[0], _HEADER)
        self.assertEqual(
            lines[1].split("\t"),
            [
                "contig-1",
                "100",
                "1800",
                "+",
                "RF01960",
                "SSU",
                "1500.0",
                "AB000001",
                "99.520",
                "1750",
                "98",
                (
                    "Eukaryota|Metazoa|Chordata|N/A|N/A|Actinopteri|N/A|"
                    "Salmoniformes|N/A|Salmonidae|N/A|Oncorhynchus|N/A|"
                    "Oncorhynchus mykiss|N/A"
                ),
                "-",
            ],
        )

    def test_crlf_input_produces_the_same_summary_as_lf_input(self) -> None:
        payload = _row().encode("utf-8")
        self.lineages.write_bytes(payload)
        self._summarize("Chordata")
        expected = self.destination.read_bytes()

        self.lineages.write_bytes(payload.replace(b"\n", b"\r\n"))
        self._summarize("Chordata")

        self.assertEqual(self.destination.read_bytes(), expected)

    def test_flags_a_hit_whose_lineage_lacks_the_host_phylum(self) -> None:
        self.lineages.write_text(
            _row(gene_id="rrna_1")
            + _row(
                gene_id="rrna_2",
                sequence="contig-2",
                subject="AB000002",
                taxid="562",
                lineage=(
                    "Bacteria;N/A;Pseudomonadota;N/A;N/A;Gammaproteobacteria;N/A;"
                    "Enterobacterales;N/A;Enterobacteriaceae;N/A;Escherichia;N/A;"
                    "Escherichia coli;N/A"
                ),
            ),
            encoding="utf-8",
        )

        lines = self._summarize("Chordata")

        self.assertEqual([line.split("\t")[-1] for line in lines[1:]], ["no", "yes"])

    def test_reports_a_hit_unresolved_at_phylum_as_non_host(self) -> None:
        self.lineages.write_text(
            _row(
                lineage=(
                    "Eukaryota;Metazoa;N/A;N/A;N/A;N/A;N/A;N/A;N/A;N/A;N/A;"
                    "N/A;N/A;N/A;N/A"
                )
            ),
            encoding="utf-8",
        )

        lines = self._summarize("Chordata")

        self.assertEqual(lines[1].split("\t")[-1], "yes")

    def test_orders_rows_by_sequence_then_coordinates(self) -> None:
        self.lineages.write_text(
            _row(gene_id="rrna_1", sequence="contig-2", start=10, end=100)
            + _row(gene_id="rrna_2", sequence="contig-1", start=900, end=2000)
            + _row(gene_id="rrna_3", sequence="contig-1", start=50, end=800),
            encoding="utf-8",
        )

        lines = self._summarize(None)

        self.assertEqual(
            [tuple(line.split("\t")[:2]) for line in lines[1:]],
            [("contig-1", "50"), ("contig-1", "900"), ("contig-2", "10")],
        )

    def test_writes_the_header_alone_when_nothing_was_classified(self) -> None:
        self.lineages.write_text("", encoding="utf-8")

        self.assertEqual(self._summarize("Chordata"), [_HEADER])

    def test_missing_host_phylum_warns_and_leaves_hits_unclassified(self) -> None:
        self.lineages.write_text(_row(), encoding="utf-8")

        lines = self._summarize(None, "3451528")

        self.assertEqual(lines[1].split("\t")[-1], "-")
        self.assertIn("3451528", self.warnings.read_text(encoding="utf-8"))

        self._summarize("Chordata", "9606")
        self.assertEqual(self.warnings.read_text(encoding="utf-8"), "")

    def test_missing_host_phylum_warns_even_without_hits(self) -> None:
        self.lineages.write_text("", encoding="utf-8")

        self.assertEqual(self._summarize(None, "3451528"), [_HEADER])
        self.assertIn("3451528", self.warnings.read_text(encoding="utf-8"))

        self._summarize(None)
        self.assertEqual(self.warnings.read_text(encoding="utf-8"), "")

    def test_rejects_a_lineage_of_unexpected_depth(self) -> None:
        self.lineages.write_text(
            _row(lineage="Eukaryota;Metazoa;Chordata"), encoding="utf-8"
        )

        with self.assertRaisesRegex(ValueError, "Invalid lineage"):
            self._summarize(None)


if __name__ == "__main__":
    unittest.main()
