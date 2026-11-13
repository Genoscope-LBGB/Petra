from __future__ import annotations

import subprocess
from collections.abc import Sequence
from pathlib import Path

_MISSING_PHYLUM = "N/A"
_MISSING_TAXID = "UNKNOWN_TAXID"
_PHYLUM_ARGUMENTS = [
    "reformat2",
    "--taxid-field",
    "1",
    "--format",
    "{phylum}",
    "--miss-rank-repl",
    _MISSING_PHYLUM,
    "--miss-taxid-repl",
    _MISSING_TAXID,
]


class HostTaxonError(RuntimeError):
    """An expected host-taxon resolution failure suitable for CLI display."""


def resolve_host_phylum(
    value: str,
    *,
    taxonomy: Path,
    taxonkit_command: Sequence[str] = ("taxonkit",),
) -> str | None:
    """Resolve a host identifier to its NCBI phylum, or None if it has none.

    A TaxID costs one taxonkit call; a scientific name costs one additional
    call to resolve the name to a TaxID. ``taxonkit_command`` may prefix
    taxonkit with a container launcher.
    """
    query = value.strip()
    taxid = (
        query
        if query.isdigit()
        else _resolve_name(query, taxonkit_command, taxonomy)
    )

    output = _run_taxonkit(
        taxonkit_command,
        _PHYLUM_ARGUMENTS,
        data_dir=taxonomy,
        input_text=f"{taxid}\n",
    )
    lines = output.splitlines()
    if len(lines) != 1:
        raise HostTaxonError(
            "taxonkit reformat2 returned an unexpected number of records for "
            f"TaxID {taxid}."
        )
    fields = lines[0].split("\t")
    if len(fields) != 2 or fields[0] != taxid:
        raise HostTaxonError(
            "taxonkit reformat2 returned an unexpected output shape for TaxID "
            f"{taxid}; expected TaxID and phylum."
        )
    phylum = fields[1].strip()
    if phylum == _MISSING_TAXID:
        raise HostTaxonError(f"TaxID {taxid} is not in the NCBI taxonomy.")
    if not phylum or phylum == _MISSING_PHYLUM:
        return None
    return phylum


def _resolve_name(
    value: str, taxonkit_command: Sequence[str], taxonomy: Path
) -> str:
    output = _run_taxonkit(
        taxonkit_command,
        ["name2taxid"],
        data_dir=taxonomy,
        input_text=f"{value}\n",
    )
    candidates: list[str] = []
    for line in output.splitlines():
        fields = line.split("\t")
        if len(fields) >= 2 and fields[1].isdigit() and fields[1] not in candidates:
            candidates.append(fields[1])
    if not candidates:
        raise HostTaxonError(
            f"Unknown host taxon {value!r}. Pass an NCBI TaxID or a scientific name."
        )
    if len(candidates) > 1:
        taxids = ", ".join(candidates)
        raise HostTaxonError(
            f"Host taxon {value!r} is ambiguous; candidate TaxIDs: {taxids}. "
            "Pass one of these TaxIDs."
        )
    return candidates[0]


def _run_taxonkit(
    command: Sequence[str],
    arguments: list[str],
    *,
    data_dir: Path,
    input_text: str,
) -> str:
    try:
        process = subprocess.run(
            [*command, "--data-dir", str(data_dir), *arguments],
            input=input_text,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as error:
        raise HostTaxonError(f"Could not run taxonkit: {error}") from error
    if process.returncode != 0:
        detail = process.stderr.strip()
        suffix = f" {detail}" if detail else ""
        raise HostTaxonError(
            f"taxonkit {arguments[0]} failed with exit code "
            f"{process.returncode}.{suffix}"
        )
    return process.stdout
