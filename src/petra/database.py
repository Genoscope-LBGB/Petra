from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path
from urllib.error import URLError
from urllib.request import Request, urlopen

NR_URL = "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
TAXON_MAP_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/"
    "prot.accession2taxid.FULL.gz"
)

_DATABASE_FILENAME = "petra.dmnd"
_NR_FILENAME = "nr.gz"
_TAXON_MAP_FILENAME = "prot.accession2taxid.FULL.gz"
_DOWNLOAD_TIMEOUT_SECONDS = 60

logger = logging.getLogger(__name__)


class DatabaseBuildError(RuntimeError):
    """An expected DIAMOND database build failure suitable for CLI display."""


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
    database_path = output_directory / _DATABASE_FILENAME
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
        raise DatabaseBuildError(f"Could not download {label} from {url}: {error}") from error
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
