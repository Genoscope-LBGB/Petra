from __future__ import annotations

import gzip
import io
import os
import stat
import subprocess
import tarfile
import sys
import tempfile
import unittest
import zipfile
from email.message import Message
from pathlib import Path
from urllib.error import HTTPError, URLError
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
    ">MH863990;k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;"
    "f__Agaricaceae;g__Agaricus;s__bisporus\nTTAACCGG\n"
)
# Rfam clan information lists family names, not accessions, and a covariance
# model declares both.
CM_FAMILIES = {
    "RF00001": "5S_rRNA",
    "RF00002": "5_8S_rRNA",
    "RF00177": "SSU_rRNA_bacteria",
    "RF01959": "SSU_rRNA_archaea",
    "RF01960": "SSU_rRNA_eukarya",
    "RF02540": "LSU_rRNA_archaea",
    "RF02541": "LSU_rRNA_bacteria",
    "RF02542": "SSU_rRNA_microsporidia",
    "RF02543": "LSU_rRNA_eukarya",
    "RF02545": "SSU_trypano_mito",
    "RF02546": "LSU_trypano_mito",
}
CM_PAYLOAD = (
    "INFERNAL1/a [1.1.4 | Dec 2020]\nNAME     {name}\nACC      {accession}\n//\n"
)
RFAM_CLANIN = (
    "CL00111\tSSU_rRNA_archaea\tSSU_rRNA_bacteria\tSSU_rRNA_eukarya\t"
    "SSU_rRNA_microsporidia\tSSU_trypano_mito\n"
    "CL00112\t5_8S_rRNA\tLSU_rRNA_archaea\tLSU_rRNA_bacteria\tLSU_rRNA_eukarya\t"
    "LSU_trypano_mito\n"
    "CL00113\t5S_rRNA\tmtPerm-5S\n"
    "CL00001\talpha_tmRNA\tbeta_tmRNA\ttmRNA\n"
)

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
        launcher.write_text(f'@"{sys.executable}" "{script}" %*\n', encoding="utf-8")
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


def _taxdump_bytes() -> bytes:
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as bundle:
        for filename in database.TAXONOMY_FILENAMES:
            payload = f"{filename}\n".encode()
            member = tarfile.TarInfo(filename)
            member.size = len(payload)
            bundle.addfile(member, io.BytesIO(payload))
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
        database.RFAM_CLANIN_URL: RFAM_CLANIN.encode("utf-8"),
    }
    for accession in database.RRNA_CM_ACCESSIONS:
        url = database.RFAM_CM_URL.format(accession=accession)
        payloads[url] = CM_PAYLOAD.format(
            name=CM_FAMILIES[accession], accession=accession
        ).encode("utf-8")
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


class InstallTaxonomyDatabaseTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.output = Path(self.temporary.name) / "database"

    def test_downloads_and_keeps_the_extracted_taxonomy(self) -> None:
        download = _fake_download({database.TAXDUMP_URL: _taxdump_bytes()})
        with patch.object(database, "_download", download):
            taxonomy = database.install_taxonomy_database(
                output=self.output,
                taxdump=None,
            )

        self.assertEqual(taxonomy, self.output / database.TAXONOMY_DIRECTORY_NAME)
        for filename in database.TAXONOMY_FILENAMES:
            self.assertEqual(
                (taxonomy / filename).read_text(encoding="utf-8"),
                f"{filename}\n",
            )
        self.assertFalse((self.output / "taxdump.tar.gz").exists())

    def test_reuses_an_existing_complete_taxonomy(self) -> None:
        taxonomy = self.output / database.TAXONOMY_DIRECTORY_NAME
        taxonomy.mkdir(parents=True)
        for filename in database.TAXONOMY_FILENAMES:
            (taxonomy / filename).touch()

        with patch.object(database, "_download") as download:
            result = database.install_taxonomy_database(
                output=self.output,
                taxdump=None,
            )

        self.assertEqual(result, taxonomy)
        download.assert_not_called()

    def test_rejects_an_incomplete_existing_taxonomy(self) -> None:
        taxonomy = self.output / database.TAXONOMY_DIRECTORY_NAME
        taxonomy.mkdir(parents=True)
        (taxonomy / "names.dmp").touch()

        with self.assertRaises(DatabaseBuildError) as raised:
            database.install_taxonomy_database(
                output=self.output,
                taxdump=None,
            )

        self.assertIn("nodes.dmp", str(raised.exception))


class _Response(io.BytesIO):
    """An HTTP response that drops the connection after ``reset_after`` bytes."""

    def __init__(self, payload: bytes, status: int, reset_after: int | None) -> None:
        super().__init__(payload)
        self.status = status
        self.reset_after = reset_after

    def read(self, size: int | None = -1) -> bytes:
        if self.reset_after is None:
            return super().read(size)
        left = self.reset_after - self.tell()
        if left <= 0:
            raise ConnectionResetError(104, "Connection reset by peer")
        if size is None or size < 0 or size > left:
            size = left
        return super().read(size)


class DownloadTest(unittest.TestCase):
    url = "https://example.org/proteins.faa.gz"
    payload = b"0123456789" * 1000

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.destination = Path(self.temporary.name) / "proteins.faa.gz"
        self.ranges: list[str | None] = []
        sleep = patch.object(database.time, "sleep")
        sleep.start()
        self.addCleanup(sleep.stop)

    def _serve(self, *outcomes: BaseException | int | None):
        """Answer each request with the next outcome: an exception to raise, a
        byte count after which the connection resets, or None to finish."""
        remaining = list(outcomes)

        def urlopen(request, timeout):
            requested = request.get_header("Range")
            self.ranges.append(requested)
            outcome = remaining.pop(0)
            if isinstance(outcome, BaseException):
                raise outcome
            offset = (
                int(requested.removeprefix("bytes=").rstrip("-")) if requested else 0
            )
            status = 206 if requested else 200
            return _Response(self.payload[offset:], status, outcome)

        return patch.object(database, "urlopen", urlopen)

    def test_retries_a_connection_reset_before_any_data(self) -> None:
        reset = URLError(ConnectionResetError(104, "Connection reset by peer"))
        with self._serve(reset, None):
            database._download(self.url, self.destination, "proteins")

        self.assertEqual(self.destination.read_bytes(), self.payload)

    def test_resumes_from_the_bytes_already_downloaded(self) -> None:
        with self._serve(3000, None):
            database._download(self.url, self.destination, "proteins")

        self.assertEqual(self.destination.read_bytes(), self.payload)
        self.assertEqual(self.ranges, [None, "bytes=3000-"])

    def test_does_not_retry_a_missing_file(self) -> None:
        missing = HTTPError(self.url, 404, "Not Found", Message(), None)
        with self._serve(missing), self.assertRaises(DatabaseBuildError) as raised:
            database._download(self.url, self.destination, "proteins")

        self.assertIn("after 1 attempt(s)", str(raised.exception))
        self.assertFalse(self.destination.exists())

    def test_gives_up_and_removes_the_partial_file(self) -> None:
        attempts = [100] * database._DOWNLOAD_ATTEMPTS
        with self._serve(*attempts), self.assertRaises(DatabaseBuildError) as raised:
            database._download(self.url, self.destination, "proteins")

        self.assertIn("Connection reset by peer", str(raised.exception))
        self.assertEqual(len(self.ranges), database._DOWNLOAD_ATTEMPTS)
        self.assertFalse(self.destination.exists())


# Captures what makedb read, since the build removes its inputs afterwards.
_DIAMOND_STUB = """\
import os
import pathlib
import shutil
import sys

arguments = sys.argv[1:]
if arguments[0] == "version":
    print("diamond version 2.2.8")
    sys.exit(0)
options = dict(zip(arguments[1::2], arguments[2::2]))
base = pathlib.Path(options["--db"])
if "--in" in options:
    shutil.copyfile(options["--in"], f"{base}.fasta")
else:
    pathlib.Path(f"{base}.fasta").write_bytes(sys.stdin.buffer.read())
shutil.copyfile(options["--taxonmap"], f"{base}.taxonmap")
exit_code = int(os.environ.get("DIAMOND_STUB_EXIT", "0"))
if exit_code == 0:
    pathlib.Path(f"{base}.dmnd").touch()
sys.exit(exit_code)
"""
NR_FASTA = ">WP_000001.1 protein [Escherichia coli]\nMKVLA\n"
NCBI_TAXON_MAP = (
    "accession\taccession.version\ttaxid\tgi\nWP_000001\tWP_000001.1\t562\t1\n"
)
MARFERRET_FASTA = ">mft0000000001 \nMSTNP\n"
MARFERRET_TAXON_MAP = (
    "accession\taccession.version\ttaxid\tgi\nNA\tmft0000000001\t2593073\tNA\n"
)


def _maybe_gunzip(payload: bytes) -> str:
    if payload[:2] == b"\x1f\x8b":
        payload = gzip.decompress(payload)
    return payload.decode("utf-8")


class BuildDiamondDatabaseMarferretTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        tools = self.directory / "tools"
        tools.mkdir()
        script = _write_stub(tools, "diamond")
        script.write_text(_DIAMOND_STUB, encoding="utf-8")
        self.output = self.directory / "database"
        patched_path = patch.dict(
            os.environ, {"PATH": f"{tools}{os.pathsep}{os.environ['PATH']}"}
        )
        patched_path.start()
        self.addCleanup(patched_path.stop)
        self.payloads = {
            database.NR_URL: _gzip_bytes(NR_FASTA),
            database.TAXON_MAP_URL: _gzip_bytes(NCBI_TAXON_MAP),
            database.MARFERRET_FASTA_URL: _gzip_bytes(MARFERRET_FASTA),
            database.MARFERRET_TAXON_MAP_URL: _gzip_bytes(MARFERRET_TAXON_MAP),
        }

    def _build(self, **inputs: Path) -> Path:
        with patch.object(database, "_download", _fake_download(self.payloads)):
            return database.build_diamond_database(
                fasta=inputs.get("fasta"),
                taxon_map=inputs.get("taxon_map"),
                output=self.output,
                marferret=True,
            )

    def _seen_by_makedb(self) -> tuple[str, str]:
        base = self.output / "petra"
        return (
            _maybe_gunzip(Path(f"{base}.fasta").read_bytes()),
            _maybe_gunzip(Path(f"{base}.taxonmap").read_bytes()),
        )

    def test_appends_marferret_to_the_downloaded_nr_inputs(self) -> None:
        self.assertTrue(self._build().is_file())

        fasta, taxon_map = self._seen_by_makedb()
        self.assertEqual(fasta, NR_FASTA + MARFERRET_FASTA)
        self.assertEqual(taxon_map, NCBI_TAXON_MAP + "NA\tmft0000000001\t2593073\tNA\n")
        leftovers = {path.name for path in self.output.iterdir()}
        self.assertEqual(leftovers, {"petra.dmnd", "petra.fasta", "petra.taxonmap"})

    def test_leaves_plain_user_inputs_untouched(self) -> None:
        fasta = self.directory / "proteins.faa"
        fasta.write_text(NR_FASTA, encoding="utf-8")
        taxon_map = self.directory / "proteins.accession2taxid"
        taxon_map.write_text(NCBI_TAXON_MAP, encoding="utf-8")

        self._build(fasta=fasta, taxon_map=taxon_map)

        seen_fasta, seen_taxon_map = self._seen_by_makedb()
        self.assertEqual(seen_fasta, NR_FASTA + MARFERRET_FASTA)
        self.assertEqual(
            seen_taxon_map, NCBI_TAXON_MAP + "NA\tmft0000000001\t2593073\tNA\n"
        )
        self.assertEqual(fasta.read_text(encoding="utf-8"), NR_FASTA)
        self.assertEqual(taxon_map.read_text(encoding="utf-8"), NCBI_TAXON_MAP)
        self.assertFalse((self.output / "petra.accession2taxid").exists())

    def test_removes_every_download_when_makedb_fails(self) -> None:
        with (
            patch.dict(os.environ, {"DIAMOND_STUB_EXIT": "1"}),
            self.assertRaises(DatabaseBuildError) as raised,
        ):
            self._build()

        self.assertIn("exit code 1", str(raised.exception))
        leftovers = {path.name for path in self.output.iterdir()}
        self.assertEqual(leftovers, {"petra.fasta", "petra.taxonmap"})


class BuildRnaDatabasesTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.tools = self.directory / "tools"
        self.tools.mkdir()
        self.output = self.directory / "database"
        self.taxonomy = self.output / database.TAXONOMY_DIRECTORY_NAME
        self.taxonomy.mkdir(parents=True)
        for filename in database.TAXONOMY_FILENAMES:
            (self.taxonomy / filename).touch()
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
            return build_rna_databases(
                output=self.output,
                taxonomy=self.taxonomy,
            )

    def test_builds_every_database_and_removes_downloads(self) -> None:
        databases = self._build()

        self.assertTrue(databases.covariance_model.is_file())
        self.assertTrue(databases.clan_info.is_file())
        self.assertTrue(databases.ssu.is_file())
        self.assertTrue(databases.lsu.is_file())
        self.assertEqual(
            databases.paths,
            (
                databases.covariance_model,
                databases.clan_info,
                databases.ssu,
                databases.lsu,
            ),
        )
        self.assertFalse((self.output / "rna_sources").exists())

    def test_uses_the_bundled_taxonomy_for_name_resolution(self) -> None:
        self._build()

        arguments = _stub_arguments(self.tools, "taxonkit")
        self.assertGreaterEqual(len(arguments), 2)
        for command in arguments:
            self.assertEqual(
                command[:3],
                ["--data-dir", str(self.taxonomy), "name2taxid"],
            )

    def test_concatenates_every_rfam_model_before_cmpress(self) -> None:
        databases = self._build()

        model = databases.covariance_model.read_text(encoding="utf-8")
        for accession, family in CM_FAMILIES.items():
            self.assertIn(f"NAME     {family}", model)
            self.assertIn(f"ACC      {accession}", model)
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
            if (
                Path(executable).stem == "makeblastdb"
                and arguments[1] == database.SSU_DATABASE_FILENAME
            ):
                index = arguments.index("-taxid_map")
                recorded["path"] = arguments[index + 1]
            return original(executable, arguments, **kwargs)

        original = database._run
        with patch.object(database, "_run", capture):
            with patch.object(database, "_download", _fake_download(_payloads())):
                with patch.object(
                    database, "_remove_source_directory", lambda path, strict: None
                ):
                    build_rna_databases(
                        output=self.output,
                        taxonomy=self.taxonomy,
                    )

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

    def test_lsu_uses_plain_headers_and_eukaryome_precedes_silva(self) -> None:
        databases = self._build()

        records = _fasta_records(databases.lsu)
        self.assertEqual(sorted(records), ["LN000001.1.2800", "MH863990"])
        self.assertEqual(records["MH863990"], "TTAACCGG")
        self.assertEqual(records["LN000001.1.2800"], "CCGGAATT")

    def test_lsu_taxid_map_combines_silva_and_resolved_eukaryome_taxids(self) -> None:
        recorded: dict[str, object] = {}

        def capture(executable: str, arguments: list[str], **kwargs: object) -> None:
            if (
                Path(executable).stem == "makeblastdb"
                and arguments[1] == database.LSU_DATABASE_FILENAME
            ):
                recorded["arguments"] = arguments
                recorded["path"] = arguments[arguments.index("-taxid_map") + 1]
            return original(executable, arguments, **kwargs)

        original = database._run
        with patch.object(database, "_run", capture):
            with patch.object(database, "_download", _fake_download(_payloads())):
                with patch.object(
                    database, "_remove_source_directory", lambda path, strict: None
                ):
                    build_rna_databases(
                        output=self.output,
                        taxonomy=self.taxonomy,
                    )

        rows = [
            line.split("\t")
            for line in Path(recorded["path"]).read_text(encoding="utf-8").splitlines()
        ]
        self.assertEqual(
            dict(rows),
            {
                "LN000001.1.2800": "9606",
                "MH863990": "4751",
            },
        )
        arguments = recorded["arguments"]
        self.assertIn("-parse_seqids", arguments)
        self.assertIn("-taxid_map", arguments)

    def test_refuses_to_overwrite_an_existing_database(self) -> None:
        for filename in (
            database.COVARIANCE_MODEL_FILENAME,
            database.CLAN_INFO_FILENAME,
            database.SSU_DATABASE_FILENAME,
            database.LSU_DATABASE_FILENAME,
        ):
            with self.subTest(filename=filename):
                self.output.mkdir(parents=True, exist_ok=True)
                (self.output / filename).touch()

                with self.assertRaises(DatabaseBuildError) as raised:
                    self._build()

                self.assertIn("RNA database already exists", str(raised.exception))
                (self.output / filename).unlink()

    def test_removes_partial_databases_when_a_tool_fails(self) -> None:
        def failing_run(
            executable: str, arguments: list[str], **kwargs: object
        ) -> None:
            if Path(executable).stem == "makeblastdb":
                raise DatabaseBuildError("makeblastdb failed with exit code 1.")

        with patch.object(database, "_run", failing_run):
            with patch.object(database, "_download", _fake_download(_payloads())):
                with self.assertRaises(DatabaseBuildError):
                    build_rna_databases(
                        output=self.output,
                        taxonomy=self.taxonomy,
                    )

        self.assertEqual(
            sorted(path.name for path in self.output.iterdir()),
            [database.TAXONOMY_DIRECTORY_NAME],
        )

    def test_reports_a_missing_executable(self) -> None:
        with patch.object(database.shutil, "which", return_value=None):
            with self.assertRaises(DatabaseBuildError) as raised:
                build_rna_databases(
                    output=self.output,
                    taxonomy=self.taxonomy,
                )

        self.assertIn("'makeblastdb' executable was not found", str(raised.exception))


class RfamClanInfoTest(unittest.TestCase):
    def test_keeps_only_multi_member_clans_of_our_families(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "Rfam.clanin"
            destination = directory / database.CLAN_INFO_FILENAME
            source.write_text(RFAM_CLANIN, encoding="utf-8")

            database._write_clan_info(
                source=source,
                destination=destination,
                families=frozenset(CM_FAMILIES.values()),
            )

            self.assertEqual(
                destination.read_text(encoding="utf-8"),
                "CL00111\tSSU_rRNA_archaea\tSSU_rRNA_bacteria\tSSU_rRNA_eukarya\t"
                "SSU_rRNA_microsporidia\tSSU_trypano_mito\n"
                "CL00112\t5_8S_rRNA\tLSU_rRNA_archaea\tLSU_rRNA_bacteria\t"
                "LSU_rRNA_eukarya\tLSU_trypano_mito\n",
            )

    def test_rejects_clan_information_listing_accessions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "Rfam.clanin"
            source.write_text(
                "CL00111\tRF00177\tRF01960\nCL00112\tRF02540\tRF02541\n",
                encoding="utf-8",
            )

            with self.assertRaises(DatabaseBuildError) as raised:
                database._write_clan_info(
                    source=source,
                    destination=directory / database.CLAN_INFO_FILENAME,
                    families=frozenset(CM_FAMILIES.values()),
                )

        self.assertIn("Rfam clan file layout changed", str(raised.exception))


class CovarianceModelNameTest(unittest.TestCase):
    def test_reads_the_family_name_declared_by_the_model(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            model = Path(temporary) / "RF01960.cm"
            model.write_text(
                CM_PAYLOAD.format(name="SSU_rRNA_eukarya", accession="RF01960"),
                encoding="utf-8",
            )

            self.assertEqual(
                database._covariance_model_name(model, "RF01960"),
                "SSU_rRNA_eukarya",
            )

    def test_rejects_a_model_of_another_accession(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            model = Path(temporary) / "RF01960.cm"
            model.write_text(
                CM_PAYLOAD.format(name="SSU_rRNA_bacteria", accession="RF00177"),
                encoding="utf-8",
            )

            with self.assertRaises(DatabaseBuildError) as raised:
                database._covariance_model_name(model, "RF01960")

        self.assertIn("does not declare a name for RF01960", str(raised.exception))

    def test_rejects_a_model_without_a_name(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            model = Path(temporary) / "RF01960.cm"
            model.write_text("INFERNAL1/a\nACC      RF01960\n//\n", encoding="utf-8")

            with self.assertRaises(DatabaseBuildError):
                database._covariance_model_name(model, "RF01960")


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
        self.assertEqual(taxon_map.accessions, {"AY855839", "JQ000001", "XX999999"})

    def test_rejects_a_taxon_map_without_accessions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "taxmap"
            path.write_text("primaryAccession\tstart\tstop\n", encoding="utf-8")

            with self.assertRaises(DatabaseBuildError) as raised:
                database._read_silva_taxon_map(path)

        self.assertIn("layout changed", str(raised.exception))


class EukaryomeRecordsTest(unittest.TestCase):
    def test_writes_plain_headers_and_looks_species_up_as_binomials(self) -> None:
        ranks, accessions, written = self._write(EUKARYOME_FASTA)

        self.assertEqual(
            ranks,
            {
                "MH863990": (
                    "Fungi",
                    "Basidiomycota",
                    "Agaricomycetes",
                    "Agaricales",
                    "Agaricaceae",
                    "Agaricus",
                    "Agaricus bisporus",
                )
            },
        )
        self.assertEqual(accessions, {"MH863990"})
        self.assertEqual(written, ">MH863990\nTTAACCGG\n")

    def test_drops_unclassified_placeholders_and_orphan_epithets(self) -> None:
        ranks, _, _ = self._write(
            ">EUK1;k__Metazoa;p__Arthropoda;c__Malacostraca;o__Amphipoda;"
            "f__Gammaridae;g__unclassified;s__unclassified\nA\n"
            ">NC_1;k___mitochondrion;p__unclassified;c__unclassified;"
            "o__unclassified;f__unclassified;g__Batrachospermum;s__unclassified\nA\n"
            ">EUK2;k__Fungi;p__unclassified;c__unclassified;o__unclassified;"
            "f__unclassified;g__unclassified;s__japonica\nA\n"
        )

        self.assertEqual(
            ranks,
            {
                "EUK1": (
                    "Metazoa",
                    "Arthropoda",
                    "Malacostraca",
                    "Amphipoda",
                    "Gammaridae",
                ),
                "NC_1": ("_mitochondrion", "Batrachospermum"),
                "EUK2": ("Fungi",),
            },
        )

    @staticmethod
    def _write(fasta: str) -> tuple[dict[str, tuple[str, ...]], set[str], str]:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "eukaryome.fasta"
            source.write_text(fasta, encoding="utf-8")
            destination = io.StringIO()
            ranks, accessions = database._write_eukaryome_records(
                source=source,
                destination=destination,
            )
        return ranks, accessions, destination.getvalue()


class ResolveTaxidsTest(unittest.TestCase):
    taxonomy = Path("database/taxonomy")

    def test_prefers_the_deepest_rank_known_to_ncbi(self) -> None:
        ranks = {"seq1": ("Eukaryota", "Obazoa", "Opisthokonta", "Unknown clade")}
        completed = subprocess.CompletedProcess(
            args=[],
            returncode=0,
            stdout="Eukaryota\t2759\nObazoa\t\nOpisthokonta\t33154\n",
            stderr="",
        )

        with patch.object(database.subprocess, "run", return_value=completed):
            taxids = database._resolve_taxids(
                ranks,
                taxonkit="taxonkit",
                taxonomy=self.taxonomy,
            )

        self.assertEqual(taxids, {"seq1": 33154})

    def test_disambiguates_a_homonym_with_the_records_accepted_ancestors(
        self,
    ) -> None:
        name2taxid = subprocess.CompletedProcess(
            args=[],
            returncode=0,
            stdout=("Eukaryota\t2759\nFungi\t4751\nDuplicata\t111\nDuplicata\t222\n"),
            stderr="",
        )
        lineages = subprocess.CompletedProcess(
            args=[],
            returncode=0,
            stdout=(
                "111\tcellular organisms;Eukaryota;Fungi;Duplicata\n"
                "222\tcellular organisms;Bacteria;Duplicata\n"
            ),
            stderr="",
        )

        with patch.object(
            database.subprocess, "run", side_effect=[name2taxid, lineages]
        ) as run:
            taxids = database._resolve_taxids(
                {"seq1": ("Eukaryota", "Fungi", "Duplicata")},
                taxonkit="taxonkit",
                taxonomy=self.taxonomy,
            )

        self.assertEqual(taxids, {"seq1": 111})
        self.assertEqual(run.call_count, 2)
        self.assertEqual(
            run.call_args_list[1].args[0],
            ["taxonkit", "--data-dir", str(self.taxonomy), "lineage"],
        )
        self.assertEqual(run.call_args_list[1].kwargs["input"], "111\n222\n")

    def test_falls_back_when_a_homonym_cannot_be_disambiguated(self) -> None:
        name2taxid = subprocess.CompletedProcess(
            args=[],
            returncode=0,
            stdout="Eukaryota\t2759\nDuplicata\t111\nDuplicata\t222\n",
            stderr="",
        )
        lineages = subprocess.CompletedProcess(
            args=[],
            returncode=0,
            stdout=(
                "111\tcellular organisms;Bacteria;Duplicata\n"
                "222\tcellular organisms;Archaea;Duplicata\n"
            ),
            stderr="",
        )

        with patch.object(
            database.subprocess, "run", side_effect=[name2taxid, lineages]
        ):
            taxids = database._resolve_taxids(
                {"seq1": ("Eukaryota", "Duplicata", "Known deeper rank")},
                taxonkit="taxonkit",
                taxonomy=self.taxonomy,
            )

        self.assertEqual(taxids, {"seq1": 2759})

    def test_name_resolution_preserves_every_distinct_candidate(self) -> None:
        completed = subprocess.CompletedProcess(
            args=[],
            returncode=0,
            stdout="Duplicata\t111\nDuplicata\t222\nDuplicata\t111\n",
            stderr="",
        )

        with patch.object(database.subprocess, "run", return_value=completed):
            candidates = database._name_to_taxid(
                ["Duplicata"],
                taxonkit="taxonkit",
                taxonomy=self.taxonomy,
            )

        self.assertEqual(candidates, {"Duplicata": [111, 222]})

    def test_skips_records_without_any_resolved_rank(self) -> None:
        ranks = {"seq1": ("Unknown clade",)}
        completed = subprocess.CompletedProcess(
            args=[], returncode=0, stdout="Unknown clade\t\n", stderr=""
        )

        with patch.object(database.subprocess, "run", return_value=completed):
            taxids = database._resolve_taxids(
                ranks,
                taxonkit="taxonkit",
                taxonomy=self.taxonomy,
            )

        self.assertEqual(taxids, {})

    def test_reports_a_taxonkit_failure(self) -> None:
        completed = subprocess.CompletedProcess(
            args=[], returncode=1, stdout="", stderr="taxonomy data not found"
        )

        with patch.object(database.subprocess, "run", return_value=completed):
            with self.assertRaises(DatabaseBuildError) as raised:
                database._resolve_taxids(
                    {"seq1": ("Fungi",)},
                    taxonkit="taxonkit",
                    taxonomy=self.taxonomy,
                )

        self.assertIn("taxonomy data not found", str(raised.exception))


if __name__ == "__main__":
    unittest.main()
