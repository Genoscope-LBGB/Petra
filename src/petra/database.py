from __future__ import annotations

import gzip
import logging
import re
import shlex
import shutil
import subprocess
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path
from urllib.error import URLError
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

NR_URL = "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
TAXON_MAP_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/"
    "prot.accession2taxid.FULL.gz"
)
RFAM_CM_URL = "https://rfam.org/family/{accession}/cm"
SILVA_SSU_URL = (
    "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/"
    "SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"
)
SILVA_SSU_TAXON_MAP_URL = (
    "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/"
    "ncbi/taxmap_embl-ebi_ena_ssu_ref_nr99_138.2.txt.gz"
)
SILVA_LSU_URL = (
    "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/"
    "SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz"
)
SILVA_LSU_TAXON_MAP_URL = (
    "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/"
    "ncbi/taxmap_embl-ebi_ena_lsu_ref_nr99_138.2.txt.gz"
)
PR2_SSU_URL = (
    "https://github.com/pr2database/pr2database/releases/download/v5.1.1/"
    "pr2_version_5.1.1_SSU_taxo_long.fasta.gz"
)
EUKARYOME_LSU_URL = (
    "https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_LSU_v2.0.zip"
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
_NR_FILENAME = "nr.gz"
_TAXON_MAP_FILENAME = "prot.accession2taxid.FULL.gz"
_CM_FILENAME = "petra.rrna.cm"
_SSU_FILENAME = "petra.ssu.fasta"
_LSU_FILENAME = "petra.lsu.fasta"
_RNA_SOURCE_DIRNAME = "rna_sources"
_CM_INDEX_SUFFIXES = (".i1f", ".i1i", ".i1m", ".i1p")
_FASTA_SUFFIXES = (".fasta", ".fa", ".fas", ".fna")
_UNIDENTIFIED_TAXID = "32644"
_PR2_COMPARTMENTS = (":apic", ":chro", ":mito", ":nucl", ":plas")
_RANK_PREFIX = re.compile(r"^[a-z]__")
_DOWNLOAD_TIMEOUT_SECONDS = 60

logger = logging.getLogger(__name__)


class DatabaseBuildError(RuntimeError):
    """An expected database build failure suitable for CLI display."""


@dataclass(frozen=True)
class RnaDatabases:
    """Databases required to detect and classify rRNA genes."""

    covariance_model: Path
    ssu: Path
    lsu: Path

    @property
    def paths(self) -> tuple[Path, Path, Path]:
        return (self.covariance_model, self.ssu, self.lsu)


@dataclass(frozen=True)
class SilvaTaxonMap:
    """SILVA taxonomy indexed by FASTA sequence identifier and by accession."""

    taxids: dict[str, int]
    unidentified: set[str]
    accessions: set[str]


def build_diamond_database(
    *,
    fasta: Path | None,
    taxon_map: Path | None,
    output: Path,
    diamond_executable: str = "diamond",
) -> Path:
    """Build ``petra.dmnd``, downloading and removing missing inputs as needed."""
    diamond = _require_diamond(diamond_executable)
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

    try:
        if fasta is None:
            _download(NR_URL, fasta_path, "NR FASTA")
            downloads.append(fasta_path)
        if taxon_map is None:
            _download(TAXON_MAP_URL, taxon_map_path, "NCBI protein taxon map")
            downloads.append(taxon_map_path)

        arguments = [
            diamond,
            "makedb",
            "--in",
            str(fasta_path),
            "--db",
            str(database_base),
            "--taxonmap",
            str(taxon_map_path),
        ]
        logger.info("Building DIAMOND database at %s", database_path)
        try:
            process = subprocess.run(arguments, cwd=output_directory, check=False)
        except OSError as error:
            raise DatabaseBuildError(f"Could not run DIAMOND: {error}") from error
        if process.returncode != 0:
            raise DatabaseBuildError(
                f"DIAMOND database build failed with exit code {process.returncode}."
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
    logger.info("Downloading %s from %s", label, url)
    logger.info("Saving %s to %s", label, destination)
    request = Request(url, headers={"User-Agent": "petra"})
    created = False
    try:
        with destination.open("xb") as output_stream:
            created = True
            with urlopen(request, timeout=_DOWNLOAD_TIMEOUT_SECONDS) as response:
                shutil.copyfileobj(response, output_stream)
    except FileExistsError as error:
        raise DatabaseBuildError(
            f"Download destination already exists: {destination}"
        ) from error
    except (OSError, URLError) as error:
        if created:
            _remove_generated_file(destination)
        raise DatabaseBuildError(
            f"Could not download {label} from {url}: {error}"
        ) from error
    logger.info("Downloaded %s to %s", label, destination)


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
    makeblastdb_executable: str = "makeblastdb",
    cmpress_executable: str = "cmpress",
    taxonkit_executable: str = "taxonkit",
) -> RnaDatabases:
    """Build the covariance model and the SSU and LSU BLAST databases."""
    makeblastdb = _require_executable(makeblastdb_executable, "BLAST+")
    cmpress = _require_executable(cmpress_executable, "Infernal")
    taxonkit = _require_executable(taxonkit_executable, "TaxonKit")

    directory = _prepare_output_directory(output)
    databases = RnaDatabases(
        covariance_model=directory / _CM_FILENAME,
        ssu=directory / _SSU_FILENAME,
        lsu=directory / _LSU_FILENAME,
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
        )
        _build_lsu_database(
            sources=sources,
            database=databases.lsu,
            makeblastdb=makeblastdb,
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
    for accession in RRNA_CM_ACCESSIONS:
        model = model_directory / f"{accession}.cm"
        _download(
            RFAM_CM_URL.format(accession=accession),
            model,
            f"Rfam covariance model {accession}",
        )
        models.append(model)

    logger.info("Concatenating %d covariance models into %s", len(models), database)
    try:
        with database.open("xb") as combined:
            for model in models:
                with model.open("rb") as stream:
                    shutil.copyfileobj(stream, combined)
    except OSError as error:
        raise DatabaseBuildError(f"Could not write {database}: {error}") from error
    _run(cmpress, [database.name], cwd=database.parent, label="cmpress")


def _build_ssu_database(
    *,
    sources: Path,
    database: Path,
    makeblastdb: str,
    taxonkit: str,
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
            for identifier, _, sequence in _iter_silva_records(
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
            for identifier, taxid in _resolve_taxids(ranks, taxonkit=taxonkit).items():
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


def _build_lsu_database(*, sources: Path, database: Path, makeblastdb: str) -> None:
    silva_fasta = sources / "silva_lsu.fasta"
    silva_taxon_map = sources / "silva_lsu.taxmap"
    eukaryome_fasta = sources / "eukaryome_lsu.fasta"
    _download_archive(SILVA_LSU_URL, silva_fasta, "SILVA LSU FASTA")
    _download_archive(SILVA_LSU_TAXON_MAP_URL, silva_taxon_map, "SILVA LSU taxon map")
    _download_archive(EUKARYOME_LSU_URL, eukaryome_fasta, "EUKARYOME LSU FASTA")

    silva = _read_silva_taxon_map(silva_taxon_map)
    logger.info("Writing %s", database)
    try:
        with database.open("x", encoding="utf-8", newline="\n") as combined:
            eukaryome_accessions = _write_eukaryome_records(
                source=eukaryome_fasta,
                destination=combined,
            )
            for identifier, lineage, sequence in _iter_silva_records(
                silva_fasta, silva.unidentified
            ):
                if _primary_accession(identifier) in eukaryome_accessions:
                    continue
                combined.write(f">{_join_lineage(identifier, lineage)}\n{sequence}\n")
    except OSError as error:
        raise DatabaseBuildError(f"Could not write {database}: {error}") from error

    _run(
        makeblastdb,
        ["-in", database.name, "-dbtype", "nucl"],
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
    return SilvaTaxonMap(taxids=taxids, unidentified=unidentified, accessions=accessions)


def _iter_silva_records(path: Path, unidentified: set[str]):
    """Yield identified SILVA records as ``(identifier, lineage, DNA sequence)``."""
    kept = 0
    for header, sequence in _iter_fasta(path):
        identifier, _, lineage = header.partition(" ")
        if not identifier or identifier in unidentified:
            continue
        kept += 1
        yield identifier, lineage, sequence.replace("U", "T")
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
    for header, sequence in _iter_fasta(source):
        for compartment in _PR2_COMPARTMENTS:
            header = header.replace(compartment, "")
        identifier, _, lineage = header.partition("|")
        identifier = identifier.strip()
        if not identifier or identifier in ranks_by_identifier:
            continue
        if _primary_accession(identifier) in excluded_accessions:
            continue
        ranks_by_identifier[identifier] = tuple(
            sys.intern(rank.strip().replace("_", " "))
            for rank in lineage.split("|")
            if rank.strip()
        )
        destination.write(f">{identifier}\n{sequence}\n")
    if not ranks_by_identifier:
        raise DatabaseBuildError(f"No usable sequences found in {source}.")
    logger.info("Kept %d PR2 sequences absent from SILVA", len(ranks_by_identifier))
    return ranks_by_identifier


def _write_eukaryome_records(*, source: Path, destination) -> set[str]:
    """Write EUKARYOME records with pipe-separated lineages, return its accessions."""
    accessions: set[str] = set()
    for header, sequence in _iter_fasta(source):
        fields = [field.strip() for field in header.split(";")]
        identifier = fields[0]
        if not identifier:
            continue
        accessions.add(_primary_accession(identifier))
        lineage = "|".join(
            rank for rank in (_RANK_PREFIX.sub("", field) for field in fields[1:]) if rank
        )
        destination.write(f">{_join_lineage(identifier, lineage)}\n{sequence}\n")
    if not accessions:
        raise DatabaseBuildError(f"No usable sequences found in {source}.")
    logger.info("Kept %d EUKARYOME accessions", len(accessions))
    return accessions


def _resolve_taxids(
    ranks_by_identifier: dict[str, tuple[str, ...]],
    *,
    taxonkit: str,
) -> dict[str, int]:
    """Map each record to the taxon of the deepest lineage rank known to NCBI."""
    names = sorted({rank for ranks in ranks_by_identifier.values() for rank in ranks})
    if not names:
        return {}
    resolved = _name_to_taxid(names, taxonkit=taxonkit)
    taxids: dict[str, int] = {}
    for identifier, ranks in ranks_by_identifier.items():
        for rank in ranks:
            taxid = resolved.get(rank)
            if taxid is not None:
                taxids[identifier] = taxid
    logger.info("Resolved a taxon for %d of %d records", len(taxids), len(ranks_by_identifier))
    return taxids


def _name_to_taxid(names: list[str], *, taxonkit: str) -> dict[str, int]:
    logger.info("Resolving %d taxon names with taxonkit", len(names))
    try:
        process = subprocess.run(
            [taxonkit, "name2taxid"],
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

    resolved: dict[str, int] = {}
    for line in process.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 2 or not fields[1]:
            continue
        try:
            taxid = int(fields[1])
        except ValueError:
            continue
        resolved.setdefault(fields[0], taxid)
    return resolved


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


def _join_lineage(identifier: str, lineage: str) -> str:
    ranks = "|".join(
        rank.replace(" ", "_")
        for rank in (part.strip() for part in lineage.replace(";", "|").split("|"))
        if rank
    )
    return f"{identifier}|{ranks}" if ranks else identifier


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
        raise DatabaseBuildError(
            f"{label} failed with exit code {process.returncode}."
        )


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
