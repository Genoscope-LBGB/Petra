from __future__ import annotations

import gzip
import io
import os
import stat
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

from petra import database
from petra.database import (
    DatabaseBuildError,
    build_rna_databases,
)

SILVA_SSU_TAXON_MAP = (
    "primaryAccession\tstart\tstop\tpath\torganism_name\ttaxid\n"
    "AY855839\t1\t1390\tBacteria;Proteobacteria;\tEscherichia coli\t562\n"
    "JQ000001\t1\t1500\tEukaryota;Fungi;\tCandida albicans\t5476\n"
    "XX999999\t1\t1200\tunidentified;\tunidentified marine\t32644\n"
)
SILVA_SSU_FASTA = (
    ">AY855839.1.1390 Bacteria;Proteobacteria;Escherichia coli\n"
    "AUGCAUGC\n"
    ">JQ000001.1.1500 Eukaryota;Fungi;Candida albicans\n"
    "GGUUAACC\n"
    ">XX999999.1.1200 unidentified;unidentified marine\n"
    "UUUUUUUU\n"
)
PR2_FASTA = (
    ">AY855839.1.1390_U|18S_rRNA|nucleus||Eukaryota|Obazoa|Opisthokonta|Fungi\n"
    "ACGTACGT\n"
    ">KX123456.1.1700_U|18S_rRNA|nucleus|strain_CBS_1|Eukaryota|Alveolata"
    "|Ciliophora|Tetrahymena_thermophila\n"
    "TTTTGGGG\n"
)
SILVA_LSU_TAXON_MAP = (
    "primaryAccession\tstart\tstop\tpath\torganism_name\ttaxid\n"
    "MH863990\t1\t3000\tEukaryota;Fungi;\tAgaricus bisporus\t5341\n"
    "LN000001\t1\t2800\tEukaryota;Metazoa;\tHomo sapiens\t9606\n"
    "ZZ000002\t1\t2600\tunidentified;\tunidentified fungus\t32644\n"
)
SILVA_LSU_FASTA = (
    ">MH863990.1.3000 Eukaryota;Fungi;Agaricus bisporus\n"
    "AAUUGGCC\n"
    ">LN000001.1.2800 Eukaryota;Metazoa;Homo sapiens\n"
    "CCGGAAUU\n"
    ">ZZ000002.1.2600 unidentified;unidentified fungus\n"
    "GGGGGGGG\n"
)
EUKARYOME_FASTA = (
    ">MH863990;k__Fungi;p__Basidiomycota;s__Agaricus_bisporus\n"
    "TTAACCGG\n"
)
CM_PAYLOAD = "INFERNAL1/a\nNAME  model\n//\n"

_STUB_TEMPLATE = """\
import pathlib
import sys

log = pathlib.Path(__file__).with_suffix(".log")
with log.open("a", encoding="utf-8") as stream:
    stream.write("\\t".join(sys.argv[1:]) + "\\n")
sys.stdout.write({stdout!r})
"""


def _write_stub(directory: Path, name: str, stdout: str = "") -> Path:
    """Install an executable stub for ``name`` that records its arguments."""
    script = directory / f"{name}.py"
    script.write_text(_STUB_TEMPLATE.format(stdout=stdout), encoding="utf-8")
    if os.name == "nt":
        launcher = directory / f"{name}.bat"
        launcher.write_text(
            f'@"{sys.executable}" "{script}" %*\n', encoding="utf-8"
        )
    else:
        launcher = directory / name
        launcher.write_text(
            f'#!/bin/sh\nexec "{sys.executable}" "{script}" "$@"\n', encoding="utf-8"
        )
        launcher.chmod(launcher.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)
    return script


def _stub_arguments(directory: Path, name: str) -> list[list[str]]:
    log = directory / f"{name}.log"
    if not log.is_file():
        return []
    return [
        line.split("\t")
        for line in log.read_text(encoding="utf-8").splitlines()
        if line
    ]


def _gzip_bytes(text: str) -> bytes:
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb") as stream:
        stream.write(text.encode("utf-8"))
    return buffer.getvalue()


def _zip_bytes(name: str, text: str) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as bundle:
        bundle.writestr(name, text)
    return buffer.getvalue()


def _payloads() -> dict[str, bytes]:
    payloads = {
        database.SILVA_SSU_URL: _gzip_bytes(SILVA_SSU_FASTA),
        database.SILVA_SSU_TAXON_MAP_URL: _gzip_bytes(SILVA_SSU_TAXON_MAP),
        database.PR2_SSU_URL: _gzip_bytes(PR2_FASTA),
        database.SILVA_LSU_URL: _gzip_bytes(SILVA_LSU_FASTA),
        database.SILVA_LSU_TAXON_MAP_URL: _gzip_bytes(SILVA_LSU_TAXON_MAP),
        database.EUKARYOME_LSU_URL: _zip_bytes(
            "General_EUK_LSU_v2.0.fasta", EUKARYOME_FASTA
        ),
    }
    for accession in database.RRNA_CM_ACCESSIONS:
        url = database.RFAM_CM_URL.format(accession=accession)
        payloads[url] = CM_PAYLOAD.replace("model", accession).encode("utf-8")
    return payloads


def _fake_download(payloads: dict[str, bytes]):
    def download(url: str, destination: Path, label: str) -> None:
        if url not in payloads:
            raise AssertionError(f"Unexpected download: {url}")
        destination.write_bytes(payloads[url])

    return download


def _fasta_records(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            identifier = line[1:]
            records[identifier] = ""
        elif identifier is not None:
            records[identifier] += line
    return records


class BuildRnaDatabasesTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.tools = self.directory / "tools"
        self.tools.mkdir()
        self.output = self.directory / "database"
        # taxonkit resolves every rank except the PR2 supergroup and division,
        # so the deepest resolved rank must win.
        _write_stub(self.tools, "makeblastdb")
        _write_stub(self.tools, "cmpress")
        _write_stub(
            self.tools,
            "taxonkit",
            stdout=(
                "Eukaryota\t2759\n"
                "Obazoa\t\n"
                "Opisthokonta\t33154\n"
                "Fungi\t4751\n"
                "Alveolata\t33630\n"
                "Ciliophora\t5878\n"
                "Tetrahymena thermophila\t312017\n"
            ),
        )
        self.patched_path = patch.dict(
            os.environ, {"PATH": f"{self.tools}{os.pathsep}{os.environ['PATH']}"}
        )
        self.patched_path.start()
        self.addCleanup(self.patched_path.stop)

    def _build(self) -> database.RnaDatabases:
        with patch.object(database, "_download", _fake_download(_payloads())):
            return build_rna_databases(output=self.output)

    def test_builds_every_database_and_removes_downloads(self) -> None:
        databases = self._build()

        self.assertTrue(databases.covariance_model.is_file())
        self.assertTrue(databases.ssu.is_file())
        self.assertTrue(databases.lsu.is_file())
        self.assertFalse((self.output / "rna_sources").exists())

    def test_concatenates_every_rfam_model_before_cmpress(self) -> None:
        databases = self._build()

        model = databases.covariance_model.read_text(encoding="utf-8")
        for accession in database.RRNA_CM_ACCESSIONS:
            self.assertIn(f"NAME  {accession}", model)
        self.assertEqual(
            _stub_arguments(self.tools, "cmpress"),
            [["petra.rrna.cm"]],
        )

    def test_ssu_drops_unidentified_silva_and_silva_duplicates_in_pr2(self) -> None:
        databases = self._build()

        records = _fasta_records(databases.ssu)
        self.assertEqual(
            sorted(records),
            ["AY855839.1.1390", "JQ000001.1.1500", "KX123456.1.1700_U"],
        )
        # SILVA sequences are stored as DNA and keep no lineage in the header.
        self.assertEqual(records["AY855839.1.1390"], "ATGCATGC")

    def test_ssu_taxid_map_uses_silva_taxids_and_deepest_pr2_rank(self) -> None:
        recorded: dict[str, str] = {}

        def capture(executable: str, arguments: list[str], **kwargs: object) -> None:
            if Path(executable).stem == "makeblastdb" and "-taxid_map" in arguments:
                index = arguments.index("-taxid_map")
                recorded["path"] = arguments[index + 1]
            return original(executable, arguments, **kwargs)

        original = database._run
        with patch.object(database, "_run", capture):
            with patch.object(database, "_download", _fake_download(_payloads())):
                with patch.object(
                    database, "_remove_source_directory", lambda path, strict: None
                ):
                    build_rna_databases(output=self.output)

        rows = [
            line.split("\t")
            for line in Path(recorded["path"]).read_text(encoding="utf-8").splitlines()
        ]
        self.assertEqual(
            dict(rows),
            {
                "AY855839.1.1390": "562",
                "JQ000001.1.1500": "5476",
                "KX123456.1.1700_U": "312017",
            },
        )

    def test_lsu_keeps_lineage_and_drops_silva_records_present_in_eukaryome(
        self,
    ) -> None:
        databases = self._build()

        headers = sorted(_fasta_records(databases.lsu))
        self.assertEqual(
            headers,
            [
                "LN000001.1.2800|Eukaryota|Metazoa|Homo_sapiens",
                "MH863990|Fungi|Basidiomycota|Agaricus_bisporus",
            ],
        )

    def test_lsu_blast_database_is_built_without_a_taxid_map(self) -> None:
        self._build()

        calls = [
            call
            for call in _stub_arguments(self.tools, "makeblastdb")
            if "petra.lsu.fasta" in call
        ]
        self.assertEqual(len(calls), 1)
        self.assertNotIn("-taxid_map", calls[0])

    def test_refuses_to_overwrite_an_existing_database(self) -> None:
        self.output.mkdir(parents=True)
        (self.output / "petra.ssu.fasta").touch()

        with self.assertRaises(DatabaseBuildError) as raised:
            self._build()

        self.assertIn("RNA database already exists", str(raised.exception))

    def test_removes_partial_databases_when_a_tool_fails(self) -> None:
        def failing_run(
            executable: str, arguments: list[str], **kwargs: object
        ) -> None:
            if Path(executable).stem == "makeblastdb":
                raise DatabaseBuildError("makeblastdb failed with exit code 1.")

        with patch.object(database, "_run", failing_run):
            with patch.object(database, "_download", _fake_download(_payloads())):
                with self.assertRaises(DatabaseBuildError):
                    build_rna_databases(output=self.output)

        self.assertEqual(sorted(path.name for path in self.output.iterdir()), [])

    def test_reports_a_missing_executable(self) -> None:
        with patch.object(database.shutil, "which", return_value=None):
            with self.assertRaises(DatabaseBuildError) as raised:
                build_rna_databases(output=self.output)

        self.assertIn("'makeblastdb' executable was not found", str(raised.exception))


class SilvaUrlTest(unittest.TestCase):
    """SILVA's 'current' symlink moved to a release without LSU exports."""

    def test_every_silva_url_pins_the_release(self) -> None:
        urls = (
            database.SILVA_SSU_URL,
            database.SILVA_SSU_TAXON_MAP_URL,
            database.SILVA_LSU_URL,
            database.SILVA_LSU_TAXON_MAP_URL,
        )
        for url in urls:
            with self.subTest(url=url):
                self.assertIn(f"release_{database.SILVA_RELEASE}/", url)
                self.assertNotIn("/current/", url)

    def test_every_silva_url_carries_the_release_in_its_filename(self) -> None:
        for url in (database.SILVA_SSU_URL, database.SILVA_LSU_URL):
            with self.subTest(url=url):
                self.assertIn(f"SILVA_{database.SILVA_RELEASE}_", url)
        for url in (
            database.SILVA_SSU_TAXON_MAP_URL,
            database.SILVA_LSU_TAXON_MAP_URL,
        ):
            with self.subTest(url=url):
                self.assertTrue(url.endswith(f"_{database.SILVA_RELEASE}.txt.gz"))


class Pr2LineageTest(unittest.TestCase):
    """A PR2 header carries gene, compartment and strain before the lineage."""

    def test_starts_at_the_domain_and_drops_the_leading_metadata(self) -> None:
        fields = [
            "18S_rRNA",
            "nucleus",
            "strain_CBS_120353",
            "Eukaryota",
            "Obazoa",
            "Knufia_epidermidis",
        ]

        self.assertEqual(
            database._pr2_lineage(fields),
            ("Eukaryota", "Obazoa", "Knufia epidermidis"),
        )

    def test_accepts_prokaryotic_domains(self) -> None:
        for domain in ("Bacteria", "Archaea"):
            with self.subTest(domain=domain):
                self.assertEqual(
                    database._pr2_lineage(["18S_rRNA", "nucleus", "", domain]),
                    (domain,),
                )

    def test_returns_nothing_without_a_domain(self) -> None:
        self.assertEqual(database._pr2_lineage(["18S_rRNA", "nucleus", ""]), ())

    def test_rejects_a_pr2_file_whose_layout_lost_the_domain(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "pr2.fasta"
            source.write_text(
                ">AB000001.1.1740_U|18S_rRNA|nucleus||Domainia|Fungi\nACGT\n",
                encoding="utf-8",
            )

            with self.assertRaises(DatabaseBuildError) as raised:
                with io.StringIO() as destination:
                    database._write_pr2_records(
                        source=source,
                        destination=destination,
                        excluded_accessions=set(),
                    )

        self.assertIn("PR2 header layout changed", str(raised.exception))


class SilvaTaxonMapTest(unittest.TestCase):
    def test_indexes_by_fasta_identifier_and_flags_unidentified(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "taxmap"
            path.write_text(SILVA_SSU_TAXON_MAP, encoding="utf-8")

            taxon_map = database._read_silva_taxon_map(path)

        self.assertEqual(taxon_map.taxids["AY855839.1.1390"], 562)
        self.assertEqual(taxon_map.unidentified, {"XX999999.1.1200"})
        self.assertEqual(
            taxon_map.accessions, {"AY855839", "JQ000001", "XX999999"}
        )

    def test_rejects_a_taxon_map_without_accessions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "taxmap"
            path.write_text("primaryAccession\tstart\tstop\n", encoding="utf-8")

            with self.assertRaises(DatabaseBuildError) as raised:
                database._read_silva_taxon_map(path)

        self.assertIn("layout changed", str(raised.exception))


class ResolveTaxidsTest(unittest.TestCase):
    def test_prefers_the_deepest_rank_known_to_ncbi(self) -> None:
        ranks = {"seq1": ("Eukaryota", "Obazoa", "Opisthokonta", "Unknown clade")}
        completed = subprocess.CompletedProcess(
            args=[],
            returncode=0,
            stdout="Eukaryota\t2759\nObazoa\t\nOpisthokonta\t33154\n",
            stderr="",
        )

        with patch.object(database.subprocess, "run", return_value=completed):
            taxids = database._resolve_taxids(ranks, taxonkit="taxonkit")

        self.assertEqual(taxids, {"seq1": 33154})

    def test_skips_records_without_any_resolved_rank(self) -> None:
        ranks = {"seq1": ("Unknown clade",)}
        completed = subprocess.CompletedProcess(
            args=[], returncode=0, stdout="Unknown clade\t\n", stderr=""
        )

        with patch.object(database.subprocess, "run", return_value=completed):
            taxids = database._resolve_taxids(ranks, taxonkit="taxonkit")

        self.assertEqual(taxids, {})

    def test_reports_a_taxonkit_failure(self) -> None:
        completed = subprocess.CompletedProcess(
            args=[], returncode=1, stdout="", stderr="taxonomy data not found"
        )

        with patch.object(database.subprocess, "run", return_value=completed):
            with self.assertRaises(DatabaseBuildError) as raised:
                database._resolve_taxids({"seq1": ("Fungi",)}, taxonkit="taxonkit")

        self.assertIn("taxonomy data not found", str(raised.exception))


if __name__ == "__main__":
    unittest.main()
