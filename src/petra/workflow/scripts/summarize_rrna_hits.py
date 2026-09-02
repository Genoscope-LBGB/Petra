import traceback
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

_LINEAGE_RANKS = 15
_LINEAGE_SEPARATOR = "|"
_NO_HOST = "-"
_NON_HOST = "yes"
_HOST = "no"
_COLUMNS = (
    "sequence_name",
    "gene_start",
    "gene_end",
    "strand",
    "rfam_accession",
    "subunit",
    "cmscan_score",
    "subject",
    "percent_identity",
    "alignment_length",
    "query_coverage",
    "lineage",
    "non_host",
)
# Field positions of a select_rrna_best_hits row extended by 'taxonkit reformat2'.
_GENE_ID = 0
_FIRST_REPORTED = 1
_GENE_START = 2
_GENE_END = 3
_LAST_REPORTED = 12
_LINEAGE = 13
_FIELDS = 14


def _non_host(ranks: list[str], host_taxon: str | None) -> str:
    """Flag a hit whose lineage does not carry the host taxon.

    The comparison is a whole-field match, so it happens at the rank the host
    taxon itself sits at. A hit resolved only above that rank carries 'N/A'
    there and is reported as non-host, which over-calls contamination rather
    than hiding it; the lineage stays in the row so it can be discounted.
    """
    if host_taxon is None:
        return _NO_HOST
    return _HOST if host_taxon in ranks else _NON_HOST


def summarize_rrna_hits(
    lineages: Path,
    destination: Path,
    *,
    host_taxon: str | None,
) -> None:
    rows: list[tuple[tuple[str, int, int, str], list[str]]] = []
    non_host_count = 0
    with lineages.open(encoding="utf-8", newline="") as input_handle:
        for line_number, line in enumerate(input_handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != _FIELDS:
                raise ValueError(
                    f"Expected {_FIELDS} fields at line {line_number} of "
                    f"{lineages}, got {len(fields)}."
                )
            ranks = fields[_LINEAGE].split(";")
            if len(ranks) != _LINEAGE_RANKS or any(not rank for rank in ranks):
                raise ValueError(
                    f"Invalid lineage at line {line_number} of {lineages}: "
                    f"expected {_LINEAGE_RANKS} non-empty semicolon-separated "
                    "ranks."
                )
            try:
                start = int(fields[_GENE_START])
                end = int(fields[_GENE_END])
            except ValueError as error:
                raise ValueError(
                    f"Invalid gene coordinates at line {line_number} of {lineages}."
                ) from error

            non_host = _non_host(ranks, host_taxon)
            if non_host == _NON_HOST:
                non_host_count += 1
            rows.append(
                (
                    (fields[_FIRST_REPORTED], start, end, fields[_GENE_ID]),
                    [
                        *fields[_FIRST_REPORTED:_LAST_REPORTED],
                        _LINEAGE_SEPARATOR.join(ranks),
                        non_host,
                    ],
                )
            )

    rows.sort(key=lambda row: row[0])
    with destination.open("w", encoding="utf-8", newline="\n") as output_handle:
        output_handle.write("\t".join(_COLUMNS) + "\n")
        for _, row in rows:
            output_handle.write("\t".join(row) + "\n")

    if not rows:
        print("No SSU or LSU gene received a BLAST hit; wrote the header only.")
    elif host_taxon is None:
        print(f"Reported {len(rows)} classified rRNA genes, no host taxon given.")
    else:
        print(
            f"Reported {len(rows)} classified rRNA genes, {non_host_count} of them "
            f"outside {host_taxon}."
        )


if "snakemake" in globals():
    workflow = globals()["snakemake"]
    log_path = Path(workflow.log[0])
    with log_path.open("w", encoding="utf-8") as log_handle:
        with redirect_stdout(log_handle), redirect_stderr(log_handle):
            try:
                summarize_rrna_hits(
                    Path(workflow.input.lineages),
                    Path(workflow.output.taxonomy),
                    host_taxon=workflow.params.host_taxon,
                )
            except Exception:
                traceback.print_exc()
                raise
