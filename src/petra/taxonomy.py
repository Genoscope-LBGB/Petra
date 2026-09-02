from __future__ import annotations

import subprocess

LINEAGE_RANKS = (
    "domain",
    "kingdom",
    "phylum",
    "subphylum",
    "superclass",
    "class",
    "superorder",
    "order",
    "suborder",
    "family",
    "subfamily",
    "genus",
    "subgenus",
    "species",
    "subspecies",
)

_LINEAGE_FORMAT = (
    "{domain|acellular root|superkingdom};{kingdom};{phylum};{subphylum};"
    "{superclass};{class};{superorder};{order};{suborder};{family};"
    "{subfamily};{genus};{subgenus};{species};{subspecies}"
)


class HostTaxonError(RuntimeError):
    """An expected host-taxon resolution failure suitable for CLI display."""


def resolve_host_taxon(
    value: str,
    *,
    taxonkit_executable: str = "taxonkit",
) -> str:
    """Resolve a user-facing taxon identifier to its canonical scientific name."""
    query = value.strip()
    taxid = query if query.isdigit() else _resolve_name(query, taxonkit_executable)

    lineage_output = _run_taxonkit(
        taxonkit_executable,
        ["lineage", "--no-lineage", "--show-name"],
        input_text=f"{taxid}\n",
    )
    lineage_lines = lineage_output.splitlines()
    if len(lineage_lines) != 1:
        raise HostTaxonError(
            "taxonkit lineage returned an unexpected number of records for "
            f"TaxID {taxid}."
        )
    lineage_fields = lineage_lines[0].split("\t")
    if len(lineage_fields) != 2 or lineage_fields[0] != taxid or not lineage_fields[1]:
        raise HostTaxonError(
            "taxonkit lineage returned an unexpected output shape for "
            f"TaxID {taxid}; expected TaxID and scientific name."
        )
    canonical_name = lineage_fields[1]

    ranked_output = _run_taxonkit(
        taxonkit_executable,
        [
            "reformat2",
            "--format",
            _LINEAGE_FORMAT,
            "--miss-rank-repl",
            "N/A",
            "--miss-taxid-repl",
            "N/A",
        ],
        input_text=f"{taxid}\n",
    )
    ranked_lines = ranked_output.splitlines()
    if len(ranked_lines) != 1:
        raise HostTaxonError(
            "taxonkit reformat2 returned an unexpected number of records for "
            f"TaxID {taxid}."
        )
    ranked_fields = ranked_lines[0].split("\t")
    if len(ranked_fields) != 2 or ranked_fields[0] != taxid:
        raise HostTaxonError(
            "taxonkit reformat2 returned an unexpected output shape for "
            f"TaxID {taxid}; expected TaxID and fixed-rank lineage."
        )
    lineage_ranks = ranked_fields[1].split(";")
    if len(lineage_ranks) != len(LINEAGE_RANKS):
        raise HostTaxonError(
            "taxonkit reformat2 returned an unexpected fixed-rank lineage for "
            f"TaxID {taxid}."
        )

    # non_host compares whole fields at exactly these ranks. Accepting a tribe or
    # another unreported rank would make every hit look non-host.
    if canonical_name not in lineage_ranks:
        ranks = ", ".join(LINEAGE_RANKS)
        raise HostTaxonError(
            f"{canonical_name} (TaxID {taxid}) sits at a rank PETRA does not "
            f"report. Choose a taxon at one of these ranks: {ranks}."
        )
    return canonical_name


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
