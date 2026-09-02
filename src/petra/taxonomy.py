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

# Ranks reformat2 can place in a lineage column of LINEAGE_RANKS: the domain
# column accepts these three spellings, every other column its own rank.
_DOMAIN_RANKS = ("domain", "acellular root", "superkingdom")
_REPORTED_RANKS = frozenset(LINEAGE_RANKS) | frozenset(_DOMAIN_RANKS)


class HostTaxonError(RuntimeError):
    """An expected host-taxon resolution failure suitable for CLI display."""


def resolve_host_taxon(
    value: str,
    *,
    taxonkit_executable: str = "taxonkit",
) -> str:
    """Resolve a user-facing taxon identifier to its canonical scientific name.

    Every taxonkit invocation reads the whole NCBI taxonomy dump, which takes
    tens of seconds on a shared filesystem, so a TaxID costs one call and a
    name two.
    """
    query = value.strip()
    taxid = query if query.isdigit() else _resolve_name(query, taxonkit_executable)

    output = _run_taxonkit(
        taxonkit_executable,
        ["lineage", "--no-lineage", "--show-name", "--show-rank"],
        input_text=f"{taxid}\n",
    )
    lines = output.splitlines()
    if len(lines) != 1:
        raise HostTaxonError(
            "taxonkit lineage returned an unexpected number of records for "
            f"TaxID {taxid}."
        )
    fields = lines[0].split("\t")
    if len(fields) != 3 or fields[0] != taxid or not fields[1]:
        raise HostTaxonError(
            "taxonkit lineage returned an unexpected output shape for TaxID "
            f"{taxid}; expected TaxID, scientific name and rank."
        )
    canonical_name, rank = fields[1], fields[2]

    # non_host compares whole lineage fields, and a lineage carries exactly the
    # ranks of LINEAGE_RANKS. A taxon at any other rank, a tribe or an unranked
    # clade, appears in no field and would make every hit look non-host.
    if rank not in _REPORTED_RANKS:
        ranks = ", ".join(LINEAGE_RANKS)
        raise HostTaxonError(
            f"{canonical_name} (TaxID {taxid}) is a {rank}, a rank PETRA does "
            f"not report. Choose a taxon at one of these ranks: {ranks}."
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
