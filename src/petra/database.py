from __future__ import annotations

import gzip
import io
import logging
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import time
import zipfile
from contextlib import ExitStack
from dataclasses import dataclass
from http.client import HTTPException
from pathlib import Path
from typing import BinaryIO
from urllib.error import HTTPError
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

NR_URL = "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
TAXON_MAP_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/"
    "prot.accession2taxid.FULL.gz"
)
RFAM_CM_URL = "https://rfam.org/family/{accession}/cm"
TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
RFAM_CLANIN_URL = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin"
# Pinned rather than tracking the 'current' symlink: SILVA 144 dropped the LSU
# exports entirely, and 138.2 is the last release providing both subunits.
SILVA_RELEASE = "138.2"
_SILVA_EXPORTS = f"https://ftp.arb-silva.de/release_{SILVA_RELEASE}/Exports"
SILVA_SSU_URL = f"{_SILVA_EXPORTS}/SILVA_{SILVA_RELEASE}_SSURef_NR99_tax_silva.fasta.gz"
SILVA_SSU_TAXON_MAP_URL = (
    f"{_SILVA_EXPORTS}/taxonomy/ncbi/"
    f"taxmap_embl-ebi_ena_ssu_ref_nr99_{SILVA_RELEASE}.txt.gz"
)
SILVA_LSU_URL = f"{_SILVA_EXPORTS}/SILVA_{SILVA_RELEASE}_LSURef_NR99_tax_silva.fasta.gz"
SILVA_LSU_TAXON_MAP_URL = (
    f"{_SILVA_EXPORTS}/taxonomy/ncbi/"
    f"taxmap_embl-ebi_ena_lsu_ref_nr99_{SILVA_RELEASE}.txt.gz"
)
PR2_SSU_URL = (
    "https://github.com/pr2database/pr2database/releases/download/v5.1.1/"
    "pr2_version_5.1.1_SSU_taxo_long.fasta.gz"
)
EUKARYOME_LSU_URL = (
    "https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_LSU_v2.0.zip"
)
MARFERRET_VERSION = "v1.1.1"
_MARFERRET_RECORD = "https://zenodo.org/records/10553848/files"
MARFERRET_FASTA_URL = (
    f"{_MARFERRET_RECORD}/MarFERReT.{MARFERRET_VERSION}.proteins.faa.gz"
)
MARFERRET_TAXON_MAP_URL = (
    f"{_MARFERRET_RECORD}/MarFERReT.{MARFERRET_VERSION}.taxonomies.tab.gz"
)

# 5S, 5.8S, and the per-domain SSU and LSU Rfam families, as used by cmscan.
RRNA_CM_ACCESSIONS = (
    "RF00001",
    "RF00002",
    "RF00177",
    "RF01959",
    "RF01960",
    "RF02540",
    "RF02541",
    "RF02542",
    "RF02543",
    "RF02545",
    "RF02546",
)

DIAMOND_DATABASE_FILENAME = "petra.dmnd"
COVARIANCE_MODEL_FILENAME = "petra.rrna.cm"
CLAN_INFO_FILENAME = "petra.rrna.clanin"
SSU_DATABASE_FILENAME = "petra.ssu.fasta"
LSU_DATABASE_FILENAME = "petra.lsu.fasta"
TAXONOMY_DIRECTORY_NAME = "taxonomy"
TAXONOMY_FILENAMES = ("names.dmp", "nodes.dmp", "delnodes.dmp", "merged.dmp")
_NR_FILENAME = "nr.gz"
_TAXON_MAP_FILENAME = "prot.accession2taxid.FULL.gz"
_MARFERRET_FASTA_FILENAME = f"MarFERReT.{MARFERRET_VERSION}.proteins.faa.gz"
_MARFERRET_TAXON_MAP_FILENAME = f"MarFERReT.{MARFERRET_VERSION}.taxonomies.tab.gz"
_COMBINED_TAXON_MAP_FILENAME = "petra.accession2taxid"
_GZIP_MAGIC = b"\x1f\x8b"
_TAXDUMP_ARCHIVE_FILENAME = "taxdump.tar.gz"
_RNA_SOURCE_DIRNAME = "rna_sources"
_CM_INDEX_SUFFIXES = (".i1f", ".i1i", ".i1m", ".i1p")
_FASTA_SUFFIXES = (".fasta", ".fa", ".fas", ".fna")
_UNIDENTIFIED_TAXID = "32644"
_PR2_COMPARTMENTS = (":apic", ":chro", ":mito", ":nucl", ":plas")
_PR2_DOMAINS = frozenset({"Eukaryota", "Bacteria", "Archaea"})
_RANK_PREFIX = re.compile(r"^[a-z]__")
_EUKARYOME_PLACEHOLDER = "unclassified"
_DOWNLOAD_TIMEOUT_SECONDS = 60
_DOWNLOAD_ATTEMPTS = 5
_DOWNLOAD_RETRY_DELAY_SECONDS = 30
# Rate limiting and gateway errors are transient; other HTTP errors are not.
_RETRYABLE_HTTP_CODES = frozenset({408, 429, 500, 502, 503, 504})

logger = logging.getLogger(__name__)


class DatabaseBuildError(RuntimeError):
    """An expected database build failure suitable for CLI display."""


@dataclass(frozen=True)
class RnaDatabases:
    """Databases required to detect and classify rRNA genes."""

    covariance_model: Path
    clan_info: Path
    ssu: Path
    lsu: Path

    @property
    def paths(self) -> tuple[Path, Path, Path, Path]:
        return (self.covariance_model, self.clan_info, self.ssu, self.lsu)


@dataclass(frozen=True)
class SilvaTaxonMap:
    """SILVA taxonomy indexed by FASTA sequence identifier and by accession."""

    taxids: dict[str, int]
    unidentified: set[str]
    accessions: set[str]


def install_taxonomy_database(*, output: Path, taxdump: Path | None) -> Path:
    """Install or reuse the NCBI taxonomy files required by TaxonKit."""
    output_directory = _prepare_output_directory(output)
    taxonomy_directory = output_directory / TAXONOMY_DIRECTORY_NAME
    if _path_exists(taxonomy_directory):
        if taxdump is not None:
            raise DatabaseBuildError(
                f"Taxonomy database already exists: {taxonomy_directory}. "
                "Omit --taxdump to reuse it."
            )
        return validate_taxonomy_directory(taxonomy_directory)

    downloaded = taxdump is None
    archive = (
        output_directory / _TAXDUMP_ARCHIVE_FILENAME
        if downloaded
        else _resolve_input(taxdump, "NCBI taxonomy dump")
    )
    if downloaded:
        _reject_existing_download(archive)

    try:
        if downloaded:
            _download(TAXDUMP_URL, archive, "NCBI taxonomy dump")
        _extract_taxonomy_database(archive, taxonomy_directory)
    except BaseException:
        _remove_source_directory(taxonomy_directory, strict=False)
        if downloaded:
            _remove_downloads([archive], strict=False)
        raise

    if downloaded:
        _remove_downloads([archive], strict=True)
    logger.info("Taxonomy database ready: %s", taxonomy_directory)
    return taxonomy_directory


def validate_taxonomy_directory(path: Path) -> Path:
    """Return a TaxonKit data directory after checking its required files."""
    missing = [
        filename for filename in TAXONOMY_FILENAMES if not (path / filename).is_file()
    ]
    if not path.is_dir() or missing:
        detail = ", ".join(missing or TAXONOMY_FILENAMES)
        raise DatabaseBuildError(f"NCBI taxonomy files not found in {path}: {detail}.")
    return path


def build_diamond_database(
    *,
    fasta: Path | None,
    taxon_map: Path | None,
    output: Path,
    taxonomy: Path,
    marferret: bool = False,
    diamond_executable: str = "diamond",
) -> Path:
    """Build ``petra.dmnd``, downloading and removing missing inputs as needed.

    ``taxonomy`` is the NCBI taxonomy directory whose ``nodes.dmp`` and
    ``names.dmp`` are embedded so searches can use ``--taxonlist`` and
    ``--taxon-exclude``. With ``marferret``, the MarFERReT proteins and taxon
    map are appended to the main inputs.
    """
    diamond = _require_diamond(diamond_executable)
    # DIAMOND runs inside the output directory, so relative paths would break.
    taxonomy = validate_taxonomy_directory(taxonomy).resolve()
    fasta_path = _resolve_input(fasta, "FASTA") if fasta is not None else None
    taxon_map_path = (
        _resolve_input(taxon_map, "taxon map") if taxon_map is not None else None
    )
    output_directory = _prepare_output_directory(output)
    database_path = output_directory / DIAMOND_DATABASE_FILENAME
    database_base = database_path.with_suffix("")

    if _path_exists(database_path):
        raise DatabaseBuildError(f"DIAMOND database already exists: {database_path}")

    downloads: list[Path] = []
    if fasta_path is None:
        fasta_path = output_directory / _NR_FILENAME
        _reject_existing_download(fasta_path)
    if taxon_map_path is None:
        taxon_map_path = output_directory / _TAXON_MAP_FILENAME
        _reject_existing_download(taxon_map_path)
    marferret_fasta = output_directory / _MARFERRET_FASTA_FILENAME
    marferret_taxon_map = output_directory / _MARFERRET_TAXON_MAP_FILENAME
    # A user-supplied map must stay untouched, so MarFERReT goes into a copy;
    # the downloaded NCBI map is ours to extend in place.
    combined_taxon_map = (
        output_directory / _COMBINED_TAXON_MAP_FILENAME
        if marferret and taxon_map is not None
        else None
    )
    if marferret:
        _reject_existing_download(marferret_fasta)
        _reject_existing_download(marferret_taxon_map)
    if combined_taxon_map is not None:
        _reject_existing_download(combined_taxon_map)

    try:
        if fasta is None:
            _download(NR_URL, fasta_path, "NR FASTA")
            downloads.append(fasta_path)
        if taxon_map is None:
            _download(TAXON_MAP_URL, taxon_map_path, "NCBI protein taxon map")
            downloads.append(taxon_map_path)
        if marferret:
            logger.info(
                "Adding MarFERReT %s to the DIAMOND database", MARFERRET_VERSION
            )
            _download(MARFERRET_FASTA_URL, marferret_fasta, "MarFERReT FASTA")
            downloads.append(marferret_fasta)
            _download(
                MARFERRET_TAXON_MAP_URL, marferret_taxon_map, "MarFERReT taxon map"
            )
            downloads.append(marferret_taxon_map)
            if combined_taxon_map is not None:
                downloads.append(combined_taxon_map)
                _copy_file(taxon_map_path, combined_taxon_map)
                taxon_map_path = combined_taxon_map
            _append_marferret_taxon_map(marferret_taxon_map, taxon_map_path)

        arguments = [diamond, "makedb"]
        if not marferret:
            arguments += ["--in", str(fasta_path)]
        arguments += [
            "--db",
            str(database_base),
            "--taxonmap",
            str(taxon_map_path),
            "--taxonnodes",
            str(taxonomy / "nodes.dmp"),
            "--taxonnames",
            str(taxonomy / "names.dmp"),
        ]
        logger.info("Building DIAMOND database at %s", database_path)
        if marferret:
            returncode = _run_makedb_on_stdin(
                arguments,
                cwd=output_directory,
                fasta=fasta_path,
                extra_fasta=marferret_fasta,
            )
        else:
            try:
                process = subprocess.run(arguments, cwd=output_directory, check=False)
            except OSError as error:
                raise DatabaseBuildError(f"Could not run DIAMOND: {error}") from error
            returncode = process.returncode
        if returncode != 0:
            raise DatabaseBuildError(
                f"DIAMOND database build failed with exit code {returncode}."
            )
        if not database_path.is_file():
            raise DatabaseBuildError(
                f"DIAMOND exited successfully but did not create {database_path}."
            )
    except BaseException:
        _remove_generated_file(database_path)
        _remove_downloads(downloads, strict=False)
        raise

    _remove_downloads(downloads, strict=True)
    logger.info("DIAMOND database ready: %s", database_path)
    return database_path


def _run_makedb_on_stdin(
    arguments: list[str], *, cwd: Path, fasta: Path, extra_fasta: Path
) -> int:
    """Stream ``fasta`` then the gzipped ``extra_fasta`` into ``diamond makedb``.

    DIAMOND takes a single ``--in`` file, and concatenating NR on disk would
    duplicate a 200 GB file. Gzip members may follow each other, so a gzipped
    FASTA is followed by the raw MarFERReT bytes; a plain one by its text.
    """
    try:
        process = subprocess.Popen(arguments, cwd=cwd, stdin=subprocess.PIPE)
    except OSError as error:
        raise DatabaseBuildError(f"Could not run DIAMOND: {error}") from error
    assert process.stdin is not None
    try:
        with process.stdin as stream:
            with fasta.open("rb") as source:
                shutil.copyfileobj(source, stream)
            _append_gzip_content(extra_fasta, stream, compressed=_is_gzip(fasta))
    except BrokenPipeError:
        # DIAMOND exited early; its exit code explains why.
        pass
    except BaseException:
        process.kill()
        process.wait()
        raise
    return process.wait()


def _is_gzip(path: Path) -> bool:
    try:
        with path.open("rb") as stream:
            return stream.read(len(_GZIP_MAGIC)) == _GZIP_MAGIC
    except OSError as error:
        raise DatabaseBuildError(f"Could not read {path}: {error}") from error


def _copy_file(source: Path, destination: Path) -> None:
    try:
        with source.open("rb") as reader, destination.open("xb") as writer:
            shutil.copyfileobj(reader, writer)
    except OSError as error:
        raise DatabaseBuildError(
            f"Could not copy {source} to {destination}: {error}"
        ) from error


def _append_gzip_content(
    source: Path, destination: BinaryIO, *, compressed: bool
) -> None:
    """Append the gzipped ``source`` to ``destination``, matching its encoding."""
    try:
        if compressed:
            with source.open("rb") as reader:
                shutil.copyfileobj(reader, destination)
            return
        with gzip.open(source, "rb") as reader:
            shutil.copyfileobj(reader, destination)
    except BrokenPipeError:
        raise
    except (OSError, EOFError) as error:
        raise DatabaseBuildError(f"Could not append {source.name}: {error}") from error


def _append_marferret_taxon_map(source: Path, destination: Path) -> None:
    """Append MarFERReT's map to ``destination`` in the destination's own columns.

    DIAMOND picks the map layout from the first header line, so every appended
    row must match it: NCBI's ``FULL`` map has two columns while MarFERReT ships
    the four-column ``accession accession.version taxid gi`` layout.
    """
    compressed = _is_gzip(destination)
    try:
        opener = gzip.open if compressed else open
        with opener(destination, "rt", encoding="utf-8", newline="") as reader:
            target_columns = reader.readline().rstrip("\r\n").split("\t")
        with ExitStack() as stack:
            reader = stack.enter_context(
                gzip.open(source, "rt", encoding="utf-8", newline="")
            )
            source_columns = reader.readline().rstrip("\r\n").split("\t")
            try:
                indices = [source_columns.index(column) for column in target_columns]
            except ValueError as error:
                raise DatabaseBuildError(
                    f"Cannot map {source.name} columns {source_columns} onto "
                    f"{destination.name} columns {target_columns}"
                ) from error
            raw = stack.enter_context(destination.open("ab"))
            writer = stack.enter_context(
                io.TextIOWrapper(
                    gzip.GzipFile(fileobj=raw, mode="wb") if compressed else raw,
                    encoding="utf-8",
                    newline="",
                )
            )
            for line in reader:
                fields = line.rstrip("\r\n").split("\t")
                writer.write("\t".join(fields[index] for index in indices) + "\n")
    except (OSError, EOFError) as error:
        raise DatabaseBuildError(f"Could not append {source.name}: {error}") from error


def _require_diamond(executable: str) -> str:
    diamond = shutil.which(executable)
    if diamond is None:
        raise DatabaseBuildError(
            "DIAMOND is required but the 'diamond' executable was not found on PATH."
        )

    try:
        process = subprocess.run(
            [diamond, "version"],
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError as error:
        raise DatabaseBuildError(f"Could not run DIAMOND: {error}") from error
    if process.returncode != 0:
        detail = process.stderr.strip() or process.stdout.strip()
        suffix = f" {detail}" if detail else ""
        raise DatabaseBuildError(
            f"DIAMOND availability check failed with exit code "
            f"{process.returncode}.{suffix}"
        )

    version = process.stdout.strip() or process.stderr.strip() or "version unknown"
    logger.info("Using DIAMOND executable %s (%s)", diamond, version)
    return diamond


def _resolve_input(path: Path, label: str) -> Path:
    try:
        resolved = path.expanduser().resolve(strict=True)
    except OSError as error:
        raise DatabaseBuildError(f"{label} file not found: {path}") from error
    if not resolved.is_file():
        raise DatabaseBuildError(f"{label} path is not a file: {resolved}")
    return resolved


def _prepare_output_directory(path: Path) -> Path:
    output = path.expanduser().resolve()
    try:
        output.mkdir(parents=True, exist_ok=True)
    except OSError as error:
        raise DatabaseBuildError(
            f"Could not create output directory {output}: {error}"
        ) from error
    if not output.is_dir():
        raise DatabaseBuildError(f"Output path is not a directory: {output}")
    return output


def _reject_existing_download(path: Path) -> None:
    if _path_exists(path):
        raise DatabaseBuildError(f"Download destination already exists: {path}")


def _path_exists(path: Path) -> bool:
    return path.exists() or path.is_symlink()


def _download(url: str, destination: Path, label: str) -> None:
    """Download ``url``, resuming after transient network failures.

    Multi-gigabyte sources (NR, MarFERReT) routinely hit connection resets, so
    each retry continues from the bytes already on disk when the server honours
    range requests.
    """
    logger.info("Downloading %s from %s", label, url)
    logger.info("Saving %s to %s", label, destination)
    try:
        destination.open("xb").close()
    except FileExistsError as error:
        raise DatabaseBuildError(
            f"Download destination already exists: {destination}"
        ) from error
    except OSError as error:
        raise DatabaseBuildError(
            f"Could not download {label} from {url}: {error}"
        ) from error

    for attempt in range(1, _DOWNLOAD_ATTEMPTS + 1):
        try:
            _download_attempt(url, destination)
        except (OSError, HTTPException) as error:
            retryable = not isinstance(error, HTTPError) or (
                error.code in _RETRYABLE_HTTP_CODES
            )
            if not retryable or attempt == _DOWNLOAD_ATTEMPTS:
                _remove_generated_file(destination)
                raise DatabaseBuildError(
                    f"Could not download {label} from {url} after {attempt} "
                    f"attempt(s): {error}"
                ) from error
            delay = _DOWNLOAD_RETRY_DELAY_SECONDS * 2 ** (attempt - 1)
            logger.warning(
                "Download of %s interrupted (%s); retrying in %d s (attempt %d/%d)",
                label,
                error,
                delay,
                attempt + 1,
                _DOWNLOAD_ATTEMPTS,
            )
            time.sleep(delay)
        except BaseException:
            _remove_generated_file(destination)
            raise
        else:
            logger.info("Downloaded %s to %s", label, destination)
            return


def _download_attempt(url: str, destination: Path) -> None:
    offset = destination.stat().st_size
    headers = {"User-Agent": "petra"}
    if offset:
        headers["Range"] = f"bytes={offset}-"
    request = Request(url, headers=headers)
    with urlopen(request, timeout=_DOWNLOAD_TIMEOUT_SECONDS) as response:
        # A server that ignores the range resends the whole file from byte 0.
        resumed = offset and response.status == 206
        with destination.open("ab" if resumed else "wb") as output_stream:
            shutil.copyfileobj(response, output_stream)


def _remove_downloads(paths: list[Path], *, strict: bool) -> None:
    failures: list[str] = []
    for path in paths:
        try:
            path.unlink()
        except FileNotFoundError:
            continue
        except OSError as error:
            logger.warning("Could not remove downloaded file %s: %s", path, error)
            failures.append(f"{path}: {error}")
        else:
            logger.info("Removed downloaded file %s", path)
    if strict and failures:
        raise DatabaseBuildError(
            "DIAMOND database was built, but downloaded input cleanup failed: "
            + "; ".join(failures)
        )


def _remove_generated_file(path: Path) -> None:
    try:
        path.unlink()
    except FileNotFoundError:
        return
    except OSError as error:
        logger.warning("Could not remove generated file %s: %s", path, error)
    else:
        logger.info("Removed incomplete generated file %s", path)


def _require_executable(executable: str, label: str) -> str:
    resolved = shutil.which(executable)
    if resolved is None:
        raise DatabaseBuildError(
            f"{label} is required but the {executable!r} executable was not found "
            "on PATH."
        )
    logger.info("Using %s executable %s", label, resolved)
    return resolved


def build_rna_databases(
    *,
    output: Path,
    taxonomy: Path,
    makeblastdb_executable: str = "makeblastdb",
    cmpress_executable: str = "cmpress",
    taxonkit_executable: str = "taxonkit",
) -> RnaDatabases:
    """Build the covariance model, clan information, and SSU/LSU databases."""
    makeblastdb = _require_executable(makeblastdb_executable, "BLAST+")
    cmpress = _require_executable(cmpress_executable, "Infernal")
    taxonkit = _require_executable(taxonkit_executable, "TaxonKit")
    taxonomy = validate_taxonomy_directory(taxonomy)

    directory = _prepare_output_directory(output)
    databases = RnaDatabases(
        covariance_model=directory / COVARIANCE_MODEL_FILENAME,
        clan_info=directory / CLAN_INFO_FILENAME,
        ssu=directory / SSU_DATABASE_FILENAME,
        lsu=directory / LSU_DATABASE_FILENAME,
    )
    for path in databases.paths:
        if _path_exists(path):
            raise DatabaseBuildError(f"RNA database already exists: {path}")
    for suffix in _CM_INDEX_SUFFIXES:
        index = _sibling(databases.covariance_model, suffix)
        if _path_exists(index):
            raise DatabaseBuildError(f"Covariance model index already exists: {index}")

    sources = directory / _RNA_SOURCE_DIRNAME
    if _path_exists(sources):
        raise DatabaseBuildError(f"Download destination already exists: {sources}")
    try:
        sources.mkdir()
    except OSError as error:
        raise DatabaseBuildError(
            f"Could not create download directory {sources}: {error}"
        ) from error

    try:
        _build_covariance_model(
            sources=sources,
            database=databases.covariance_model,
            cmpress=cmpress,
        )
        _build_ssu_database(
            sources=sources,
            database=databases.ssu,
            makeblastdb=makeblastdb,
            taxonkit=taxonkit,
            taxonomy=taxonomy,
        )
        _build_lsu_database(
            sources=sources,
            database=databases.lsu,
            makeblastdb=makeblastdb,
            taxonkit=taxonkit,
            taxonomy=taxonomy,
        )
    except BaseException:
        for path in databases.paths:
            _remove_database_artifacts(path)
        _remove_source_directory(sources, strict=False)
        raise

    _remove_source_directory(sources, strict=True)
    for path in databases.paths:
        logger.info("RNA database ready: %s", path)
    return databases


def _build_covariance_model(*, sources: Path, database: Path, cmpress: str) -> None:
    model_directory = sources / "rfam"
    try:
        model_directory.mkdir()
    except OSError as error:
        raise DatabaseBuildError(
            f"Could not create download directory {model_directory}: {error}"
        ) from error

    models = []
    families = []
    for accession in RRNA_CM_ACCESSIONS:
        model = model_directory / f"{accession}.cm"
        _download(
            RFAM_CM_URL.format(accession=accession),
            model,
            f"Rfam covariance model {accession}",
        )
        models.append(model)
        families.append(_covariance_model_name(model, accession))

    logger.info("Concatenating %d covariance models into %s", len(models), database)
    try:
        with database.open("xb") as combined:
            for model in models:
                with model.open("rb") as stream:
                    shutil.copyfileobj(stream, combined)
    except OSError as error:
        raise DatabaseBuildError(f"Could not write {database}: {error}") from error
    _run(cmpress, [database.name], cwd=database.parent, label="cmpress")
    clan_source = sources / "Rfam.clanin"
    _download(RFAM_CLANIN_URL, clan_source, "Rfam clan information")
    # cmscan uses --clanin to mark overlapping SSU/LSU models as competitors;
    # filtering prevents references to models absent from our combined CM.
    _write_clan_info(
        source=clan_source,
        destination=database.with_name(CLAN_INFO_FILENAME),
        families=frozenset(families),
    )


def _covariance_model_name(model: Path, accession: str) -> str:
    """Read the family name a covariance model is known by.

    Rfam clan information lists family names rather than accessions, so the
    names are taken from the downloaded models instead of from a second list
    that would have to be kept in step with RRNA_CM_ACCESSIONS.
    """
    name: str | None = None
    declared_accession: str | None = None
    try:
        with model.open(encoding="utf-8", errors="replace") as stream:
            for line in stream:
                key, _, value = line.partition(" ")
                if key == "NAME" and name is None:
                    name = value.strip()
                elif key == "ACC" and declared_accession is None:
                    declared_accession = value.strip()
                if name and declared_accession:
                    break
    except OSError as error:
        raise DatabaseBuildError(f"Could not read {model}: {error}") from error
    if not name or declared_accession != accession:
        raise DatabaseBuildError(
            f"{model} does not declare a name for {accession}; the Rfam "
            "covariance model layout changed."
        )
    return name


def _write_clan_info(
    *,
    source: Path,
    destination: Path,
    families: frozenset[str],
) -> None:
    clans: list[tuple[str, list[str]]] = []
    try:
        with source.open(encoding="utf-8", errors="replace") as stream:
            for line in stream:
                fields = line.split()
                if len(fields) < 3:
                    continue
                members = list(
                    dict.fromkeys(family for family in fields[1:] if family in families)
                )
                if len(members) >= 2:
                    clans.append((fields[0], members))
    except OSError as error:
        raise DatabaseBuildError(f"Could not read {source}: {error}") from error

    if not clans:
        raise DatabaseBuildError(
            f"No clan of the covariance model families found in {source}; the "
            "Rfam clan file layout changed."
        )

    try:
        with destination.open("x", encoding="utf-8", newline="\n") as stream:
            for clan, members in clans:
                stream.write("\t".join((clan, *members)) + "\n")
    except OSError as error:
        raise DatabaseBuildError(f"Could not write {destination}: {error}") from error
    logger.info(
        "Kept %d Rfam clans containing %d families",
        len(clans),
        sum(len(members) for _, members in clans),
    )


def _build_ssu_database(
    *,
    sources: Path,
    database: Path,
    makeblastdb: str,
    taxonkit: str,
    taxonomy: Path,
) -> None:
    silva_fasta = sources / "silva_ssu.fasta"
    silva_taxon_map = sources / "silva_ssu.taxmap"
    pr2_fasta = sources / "pr2_ssu.fasta"
    _download_archive(SILVA_SSU_URL, silva_fasta, "SILVA SSU FASTA")
    _download_archive(SILVA_SSU_TAXON_MAP_URL, silva_taxon_map, "SILVA SSU taxon map")
    _download_archive(PR2_SSU_URL, pr2_fasta, "PR2 SSU FASTA")

    silva = _read_silva_taxon_map(silva_taxon_map)
    taxid_map = sources / "ssu.taxid_map"
    logger.info("Writing %s", database)
    try:
        with (
            database.open("x", encoding="utf-8", newline="\n") as combined,
            taxid_map.open("x", encoding="utf-8", newline="\n") as taxids,
        ):
            for identifier, sequence in _iter_silva_records(
                silva_fasta, silva.unidentified
            ):
                combined.write(f">{identifier}\n{sequence}\n")
                taxid = silva.taxids.get(identifier)
                if taxid is not None:
                    taxids.write(f"{identifier}\t{taxid}\n")
            ranks = _write_pr2_records(
                source=pr2_fasta,
                destination=combined,
                excluded_accessions=silva.accessions,
            )
            for identifier, taxid in _resolve_taxids(
                ranks, taxonkit=taxonkit, taxonomy=taxonomy
            ).items():
                taxids.write(f"{identifier}\t{taxid}\n")
    except OSError as error:
        raise DatabaseBuildError(f"Could not write {database}: {error}") from error

    _run(
        makeblastdb,
        [
            "-in",
            database.name,
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-taxid_map",
            str(taxid_map),
        ],
        cwd=database.parent,
        label="makeblastdb",
    )


def _build_lsu_database(
    *,
    sources: Path,
    database: Path,
    makeblastdb: str,
    taxonkit: str,
    taxonomy: Path,
) -> None:
    silva_fasta = sources / "silva_lsu.fasta"
    silva_taxon_map = sources / "silva_lsu.taxmap"
    eukaryome_fasta = sources / "eukaryome_lsu.fasta"
    _download_archive(SILVA_LSU_URL, silva_fasta, "SILVA LSU FASTA")
    _download_archive(SILVA_LSU_TAXON_MAP_URL, silva_taxon_map, "SILVA LSU taxon map")
    _download_archive(EUKARYOME_LSU_URL, eukaryome_fasta, "EUKARYOME LSU FASTA")

    silva = _read_silva_taxon_map(silva_taxon_map)
    taxid_map = sources / "lsu.taxid_map"
    logger.info("Writing %s", database)
    try:
        with (
            database.open("x", encoding="utf-8", newline="\n") as combined,
            taxid_map.open("x", encoding="utf-8", newline="\n") as taxids,
        ):
            ranks, eukaryome_accessions = _write_eukaryome_records(
                source=eukaryome_fasta,
                destination=combined,
            )
            for identifier, sequence in _iter_silva_records(
                silva_fasta, silva.unidentified
            ):
                if _primary_accession(identifier) in eukaryome_accessions:
                    continue
                combined.write(f">{identifier}\n{sequence}\n")
                taxid = silva.taxids.get(identifier)
                if taxid is not None:
                    taxids.write(f"{identifier}\t{taxid}\n")
            resolved = _resolve_taxids(ranks, taxonkit=taxonkit, taxonomy=taxonomy)
            for identifier, taxid in resolved.items():
                taxids.write(f"{identifier}\t{taxid}\n")
    except OSError as error:
        raise DatabaseBuildError(f"Could not write {database}: {error}") from error

    _run(
        makeblastdb,
        [
            "-in",
            database.name,
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-taxid_map",
            str(taxid_map),
        ],
        cwd=database.parent,
        label="makeblastdb",
    )


def _read_silva_taxon_map(path: Path) -> SilvaTaxonMap:
    """Index a SILVA taxonomy map by the sequence identifiers used in its FASTA."""
    taxids: dict[str, int] = {}
    unidentified: set[str] = set()
    accessions: set[str] = set()
    try:
        with path.open(encoding="utf-8", errors="replace") as stream:
            stream.readline()
            for line in stream:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6:
                    continue
                accession, start, stop = fields[0], fields[1], fields[2]
                if not accession:
                    continue
                accessions.add(accession)
                identifier = f"{accession}.{start}.{stop}"
                if fields[5] == _UNIDENTIFIED_TAXID:
                    unidentified.add(identifier)
                    continue
                try:
                    taxids[identifier] = int(fields[5])
                except ValueError:
                    continue
    except OSError as error:
        raise DatabaseBuildError(f"Could not read {path}: {error}") from error
    if not accessions:
        raise DatabaseBuildError(
            f"No accessions found in {path}; the SILVA taxonomy map layout changed."
        )
    logger.info(
        "Read %d SILVA accessions from %s, %d of them unidentified",
        len(accessions),
        path,
        len(unidentified),
    )
    return SilvaTaxonMap(
        taxids=taxids, unidentified=unidentified, accessions=accessions
    )


def _iter_silva_records(path: Path, unidentified: set[str]):
    """Yield identified SILVA records as ``(identifier, DNA sequence)``."""
    kept = 0
    for header, sequence in _iter_fasta(path):
        identifier = header.partition(" ")[0]
        if not identifier or identifier in unidentified:
            continue
        kept += 1
        yield identifier, sequence.replace("U", "T")
    if not kept:
        raise DatabaseBuildError(f"No identified sequences found in {path}.")


def _write_pr2_records(
    *,
    source: Path,
    destination,
    excluded_accessions: set[str],
) -> dict[str, tuple[str, ...]]:
    """Write PR2 records absent from SILVA and return their lineages per record."""
    ranks_by_identifier: dict[str, tuple[str, ...]] = {}
    without_domain = 0
    for header, sequence in _iter_fasta(source):
        for compartment in _PR2_COMPARTMENTS:
            header = header.replace(compartment, "")
        fields = [field.strip() for field in header.split("|")]
        identifier = fields[0]
        if not identifier or identifier in ranks_by_identifier:
            continue
        if _primary_accession(identifier) in excluded_accessions:
            continue
        ranks = _pr2_lineage(fields[1:])
        if not ranks:
            without_domain += 1
            continue
        ranks_by_identifier[identifier] = ranks
        destination.write(f">{identifier}\n{sequence}\n")
    if not ranks_by_identifier:
        raise DatabaseBuildError(
            f"No sequence in {source} carries a recognizable domain; the PR2 "
            "header layout changed."
        )
    if without_domain:
        logger.warning(
            "Skipped %d PR2 sequences without a recognizable domain", without_domain
        )
    logger.info("Kept %d PR2 sequences absent from SILVA", len(ranks_by_identifier))
    return ranks_by_identifier


def _pr2_lineage(fields: list[str]) -> tuple[str, ...]:
    """Return the lineage of a PR2 header, from its domain field onwards.

    A PR2 header interleaves non-taxonomic fields before the lineage, currently
    gene, compartment and strain. Anchoring on the domain keeps the gene name
    and the strain label out of the taxon names handed to taxonkit.
    """
    for position, field in enumerate(fields):
        if field in _PR2_DOMAINS:
            return tuple(
                sys.intern(rank.replace("_", " ")) for rank in fields[position:] if rank
            )
    return ()


def _write_eukaryome_records(
    *,
    source: Path,
    destination,
) -> tuple[dict[str, tuple[str, ...]], set[str]]:
    """Write EUKARYOME records and return their searchable ranks and accessions."""
    ranks_by_identifier: dict[str, tuple[str, ...]] = {}
    accessions: set[str] = set()
    for header, sequence in _iter_fasta(source):
        fields = [field.strip() for field in header.split(";")]
        identifier = fields[0]
        if not identifier or identifier in ranks_by_identifier:
            continue
        ranks_by_identifier[identifier] = _eukaryome_ranks(fields[1:])
        accessions.add(_primary_accession(identifier))
        destination.write(f">{identifier}\n{sequence}\n")
    if not ranks_by_identifier:
        raise DatabaseBuildError(f"No usable sequences found in {source}.")
    logger.info("Kept %d EUKARYOME accessions", len(accessions))
    return ranks_by_identifier, accessions


def _eukaryome_ranks(fields: list[str]) -> tuple[str, ...]:
    """Return the names of a EUKARYOME lineage that are safe to look up in NCBI.

    TaxonKit matches names case-insensitively and against every name class, so
    two EUKARYOME conventions resolve to unrelated taxa: the 'unclassified'
    filler matches NCBI's 'unclassified sequences' (no ranked lineage), and the
    bare species epithet matches homonymous genera ('japonica' -> Japonica).
    Species are therefore looked up as binomials, and only under a named genus.
    """
    ranks: list[str] = []
    genus = ""
    for field in fields:
        name = _RANK_PREFIX.sub("", field)
        if not name or name == _EUKARYOME_PLACEHOLDER:
            continue
        if field.startswith("g__"):
            genus = name
        elif field.startswith("s__"):
            if not genus:
                continue
            name = f"{genus} {name}"
        ranks.append(name)
    return tuple(ranks)


def _resolve_taxids(
    ranks_by_identifier: dict[str, tuple[str, ...]],
    *,
    taxonkit: str,
    taxonomy: Path,
) -> dict[str, int]:
    """Map each record to the deepest unambiguous NCBI taxon in its lineage."""
    names = sorted({rank for ranks in ranks_by_identifier.values() for rank in ranks})
    if not names:
        return {}
    candidates_by_name = _name_to_taxid(names, taxonkit=taxonkit, taxonomy=taxonomy)
    ambiguous_taxids = {
        taxid
        for candidates in candidates_by_name.values()
        if len(candidates) > 1
        for taxid in candidates
    }
    lineages = _taxid_lineage_names(
        ambiguous_taxids, taxonkit=taxonkit, taxonomy=taxonomy
    )

    taxids: dict[str, int] = {}
    for identifier, ranks in ranks_by_identifier.items():
        accepted_names: list[str] = []
        for rank in ranks:
            candidates = candidates_by_name.get(rank, [])
            if not candidates:
                continue
            if len(candidates) == 1:
                taxid = candidates[0]
            else:
                matching = [
                    candidate
                    for candidate in candidates
                    if all(
                        ancestor in lineages.get(candidate, set())
                        for ancestor in accepted_names
                    )
                ]
                if len(matching) != 1:
                    break
                taxid = matching[0]
            taxids[identifier] = taxid
            accepted_names.append(rank)
    logger.info(
        "Saw %d ambiguous taxon names; resolved a taxon for %d of %d records",
        sum(len(candidates) > 1 for candidates in candidates_by_name.values()),
        len(taxids),
        len(ranks_by_identifier),
    )
    return taxids


def _name_to_taxid(
    names: list[str], *, taxonkit: str, taxonomy: Path
) -> dict[str, list[int]]:
    logger.info("Resolving %d taxon names with taxonkit", len(names))
    try:
        process = subprocess.run(
            [taxonkit, "--data-dir", str(taxonomy), "name2taxid"],
            input="\n".join(names) + "\n",
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as error:
        raise DatabaseBuildError(f"Could not run taxonkit: {error}") from error
    if process.returncode != 0:
        detail = process.stderr.strip()
        suffix = f" {detail}" if detail else ""
        raise DatabaseBuildError(
            f"taxonkit name2taxid failed with exit code {process.returncode}.{suffix}"
        )

    resolved: dict[str, list[int]] = {}
    for line in process.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 2 or not fields[1]:
            continue
        try:
            taxid = int(fields[1])
        except ValueError:
            continue
        candidates = resolved.setdefault(fields[0], [])
        if taxid not in candidates:
            candidates.append(taxid)
    return resolved


def _taxid_lineage_names(
    taxids: set[int], *, taxonkit: str, taxonomy: Path
) -> dict[int, set[str]]:
    if not taxids:
        return {}
    ordered_taxids = sorted(taxids)
    logger.info("Reading lineages for %d ambiguous taxon candidates", len(taxids))
    try:
        process = subprocess.run(
            [taxonkit, "--data-dir", str(taxonomy), "lineage"],
            input="\n".join(str(taxid) for taxid in ordered_taxids) + "\n",
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as error:
        raise DatabaseBuildError(f"Could not run taxonkit: {error}") from error
    if process.returncode != 0:
        detail = process.stderr.strip()
        suffix = f" {detail}" if detail else ""
        raise DatabaseBuildError(
            f"taxonkit lineage failed with exit code {process.returncode}.{suffix}"
        )

    lineages: dict[int, set[str]] = {}
    for line in process.stdout.splitlines():
        fields = line.split("\t")
        try:
            taxid = int(fields[0])
        except (IndexError, ValueError) as error:
            raise DatabaseBuildError(
                "Could not read taxonkit lineage output; its layout changed."
            ) from error
        if len(fields) < 2:
            raise DatabaseBuildError(
                "Could not read taxonkit lineage output; its layout changed."
            )
        lineages[taxid] = {
            name.strip() for name in fields[1].split(";") if name.strip()
        }
    if not lineages:
        raise DatabaseBuildError(
            "Could not read taxonkit lineage output; its layout changed."
        )
    return lineages


def _iter_fasta(path: Path):
    try:
        with path.open(encoding="utf-8", errors="replace") as stream:
            header: str | None = None
            chunks: list[str] = []
            for line in stream:
                line = line.strip()
                if line.startswith(">"):
                    if header is not None:
                        yield header, "".join(chunks)
                    header = line[1:]
                    chunks = []
                elif line and header is not None:
                    chunks.append(line)
            if header is not None:
                yield header, "".join(chunks)
    except OSError as error:
        raise DatabaseBuildError(f"Could not read {path}: {error}") from error


def _primary_accession(identifier: str) -> str:
    return identifier.split(".", 1)[0]


def _download_archive(url: str, destination: Path, label: str) -> None:
    suffix = Path(urlsplit(url).path).suffix.lower()
    archive = _sibling(destination, suffix) if suffix else destination
    _download(url, archive, label)
    if archive == destination:
        return
    try:
        if suffix == ".zip":
            _extract_fasta_member(archive, destination)
        else:
            with (
                gzip.open(archive, "rb") as compressed,
                destination.open("xb") as extracted,
            ):
                shutil.copyfileobj(compressed, extracted)
    except (OSError, zipfile.BadZipFile) as error:
        _remove_generated_file(destination)
        _remove_generated_file(archive)
        raise DatabaseBuildError(
            f"Could not extract {label} from {archive}: {error}"
        ) from error
    except BaseException:
        _remove_generated_file(destination)
        _remove_generated_file(archive)
        raise
    logger.info("Extracted %s to %s", label, destination)
    _remove_downloads([archive], strict=False)


def _extract_taxonomy_database(archive: Path, destination: Path) -> None:
    logger.info("Extracting NCBI taxonomy dump from %s", archive)
    try:
        with tarfile.open(archive, mode="r:*") as bundle:
            members_by_name: dict[str, tarfile.TarInfo] = {}
            for member in bundle.getmembers():
                name = Path(member.name).name
                if name not in TAXONOMY_FILENAMES or not member.isfile():
                    continue
                if name in members_by_name:
                    raise DatabaseBuildError(
                        f"NCBI taxonomy dump contains more than one {name} file: "
                        f"{archive}"
                    )
                members_by_name[name] = member

            missing = [
                filename
                for filename in TAXONOMY_FILENAMES
                if filename not in members_by_name
            ]
            if missing:
                raise DatabaseBuildError(
                    f"NCBI taxonomy dump is missing required files: "
                    f"{', '.join(missing)}."
                )

            destination.mkdir()
            for filename in TAXONOMY_FILENAMES:
                source = bundle.extractfile(members_by_name[filename])
                if source is None:
                    raise DatabaseBuildError(
                        f"Could not read {filename} from NCBI taxonomy dump."
                    )
                with source, (destination / filename).open("xb") as extracted:
                    shutil.copyfileobj(source, extracted)
    except DatabaseBuildError:
        raise
    except (OSError, tarfile.TarError) as error:
        raise DatabaseBuildError(
            f"Could not extract NCBI taxonomy dump {archive}: {error}"
        ) from error


def _extract_fasta_member(archive: Path, destination: Path) -> None:
    with zipfile.ZipFile(archive) as bundle:
        members = [
            member
            for member in bundle.infolist()
            if not member.is_dir()
            and Path(member.filename).suffix.lower() in _FASTA_SUFFIXES
        ]
        if len(members) != 1:
            expected = ", ".join(_FASTA_SUFFIXES)
            raise DatabaseBuildError(
                f"Expected exactly one file with a {expected} extension in "
                f"{archive}, found {len(members)}."
            )
        with (
            bundle.open(members[0]) as member,
            destination.open("xb") as extracted,
        ):
            shutil.copyfileobj(member, extracted)


def _run(executable: str, arguments: list[str], *, cwd: Path, label: str) -> None:
    command = [executable, *arguments]
    logger.info("Running %s", shlex.join(command))
    try:
        process = subprocess.run(command, cwd=cwd, check=False)
    except OSError as error:
        raise DatabaseBuildError(f"Could not run {label}: {error}") from error
    if process.returncode != 0:
        raise DatabaseBuildError(f"{label} failed with exit code {process.returncode}.")


def _sibling(path: Path, suffix: str) -> Path:
    return path.with_name(path.name + suffix)


def _remove_database_artifacts(path: Path) -> None:
    _remove_generated_file(path)
    for artifact in sorted(path.parent.glob(f"{path.name}.*")):
        _remove_generated_file(artifact)


def _remove_source_directory(path: Path, *, strict: bool) -> None:
    shutil.rmtree(path, ignore_errors=True)
    if path.exists():
        logger.warning("Could not remove download directory %s", path)
        if strict:
            raise DatabaseBuildError(
                "RNA databases were built, but downloaded input cleanup failed: "
                f"{path} could not be removed."
            )
        return
    logger.info("Removed download directory %s", path)
