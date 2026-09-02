from __future__ import annotations

import subprocess

_MISSING_PHYLUM = "N/A"
_PHYLUM_ARGUMENTS = [
    "reformat2",
    "--taxid-field",
    "1",
    "--format",
    "{phylum}",
    "--miss-rank-repl",
    _MISSING_PHYLUM,
    "--miss-taxid-repl",
    _MISSING_PHYLUM,
]


class HostTaxonError(RuntimeError):
    """An expected host-taxon resolution failure suitable for CLI display."""


def resolve_host_phylum(
    value: str,
    *,
    taxonkit_executable: str = "taxonkit",
) -> str:
    """Resolve a user-facing host identifier to its NCBI phylum name.

    A TaxID costs one taxonkit call; a scientific name costs one additional
    call to resolve the name to a TaxID.
    """
    query = value.strip()
    taxid = query if query.isdigit() else _resolve_name(query, taxonkit_executable)

    output = _run_taxonkit(
        taxonkit_executable,
        _PHYLUM_ARGUMENTS,
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
    if not phylum or phylum == _MISSING_PHYLUM:
        raise HostTaxonError(f"TaxID {taxid} has no phylum in the NCBI taxonomy.")
    return phylum


def _resolve_name(value: str, taxonkit_executable: str) -> str:
    output = _run_taxonkit(
        taxonkit_executable,
        ["name2taxid"],
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
    executable: str,
    arguments: list[str],
    *,
    input_text: str,
) -> str:
    try:
        process = subprocess.run(
            [executable, *arguments],
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
