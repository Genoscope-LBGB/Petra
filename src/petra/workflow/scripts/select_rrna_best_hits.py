import csv
import traceback
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass
from pathlib import Path

# Field positions of a 'blastn -outfmt "6 std qlen slen qcovs staxid"' row.
_QUERY = 0
_SUBJECT = 1
_IDENTITY = 2
_ALIGNMENT_LENGTH = 3
_BITSCORE = 11
_QUERY_COVERAGE = 14
_TAXID = 15
_MINIMUM_FIELDS = 16
_MISSING_TAXID = "0"

_GENE_COLUMNS = (
    "gene_id",
    "sequence_name",
    "gene_start",
    "gene_end",
    "strand",
    "rfam_accession",
    "subunit",
    "cmscan_score",
)


@dataclass(frozen=True)
class BlastHit:
    """One BLAST hit, keeping the reported numbers as printed."""

    subject: str
    identity: str
    alignment_length: str
    query_coverage: str
    taxid: str
    bitscore_value: float
    identity_value: float
    alignment_length_value: int

    @property
    def rank(self) -> tuple[float, float, int, str]:
        """Sort key electing the best hit; ties break deterministically.

        BLAST emits equally scoring subjects in an order that is not stable
        across runs, so identity, alignment length and finally the subject
        accession decide instead of whichever line came first.
        """
        return (
            -self.bitscore_value,
            -self.identity_value,
            -self.alignment_length_value,
            self.subject,
        )


def _taxid(value: str) -> str:
    """Normalize the staxid field, which is empty for subjects without a TaxID."""
    candidate = value.split(";", 1)[0].strip()
    return candidate if candidate.isdecimal() else _MISSING_TAXID


def _best_hits(source: Path, *, minimum_alignment_length: int) -> dict[str, BlastHit]:
    best: dict[str, BlastHit] = {}
    with source.open(encoding="utf-8", newline="") as input_handle:
        for line_number, line in enumerate(input_handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < _MINIMUM_FIELDS:
                raise ValueError(
                    f"Expected at least {_MINIMUM_FIELDS} fields at line "
                    f"{line_number} of {source}, got {len(fields)}. blastn must "
                    "run with '-outfmt \"6 std qlen slen qcovs staxid\"'."
                )
            try:
                bitscore = float(fields[_BITSCORE])
                identity = float(fields[_IDENTITY])
                alignment_length = int(fields[_ALIGNMENT_LENGTH])
            except ValueError as error:
                raise ValueError(
                    f"Invalid alignment metrics at line {line_number} of {source}."
                ) from error
            if alignment_length < minimum_alignment_length:
                continue
            hit = BlastHit(
                subject=fields[_SUBJECT],
                identity=fields[_IDENTITY],
                alignment_length=fields[_ALIGNMENT_LENGTH],
                query_coverage=fields[_QUERY_COVERAGE],
                taxid=_taxid(fields[_TAXID]),
                bitscore_value=bitscore,
                identity_value=identity,
                alignment_length_value=alignment_length,
            )
            query = fields[_QUERY]
            current = best.get(query)
            if current is None or hit.rank < current.rank:
                best[query] = hit
    return best


def _genes(source: Path) -> list[dict[str, str]]:
    with source.open(encoding="utf-8", newline="") as input_handle:
        reader = csv.reader(input_handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError(f"Empty gene table {source}.") from error
        if tuple(header) != _GENE_COLUMNS:
            raise ValueError(
                f"Unexpected gene table columns in {source}: {header}. Expected "
                f"{list(_GENE_COLUMNS)}."
            )
        return [dict(zip(_GENE_COLUMNS, row, strict=True)) for row in reader if row]


def select_rrna_best_hits(
    genes_table: Path,
    ssu_alignment: Path,
    lsu_alignment: Path,
    destination: Path,
    *,
    minimum_alignment_length: int,
) -> None:
    if minimum_alignment_length < 1:
        raise ValueError("The minimum alignment length must be at least 1.")

    genes = _genes(genes_table)
    hits = {
        "SSU": _best_hits(
            ssu_alignment, minimum_alignment_length=minimum_alignment_length
        ),
        "LSU": _best_hits(
            lsu_alignment, minimum_alignment_length=minimum_alignment_length
        ),
    }

    retained = 0
    with destination.open("w", encoding="utf-8", newline="\n") as output_handle:
        for gene in genes:
            subunit = gene["subunit"]
            if subunit not in hits:
                raise ValueError(
                    f"Gene {gene['gene_id']} carries an unknown subunit "
                    f"{subunit!r} in {genes_table}."
                )
            hit = hits[subunit].get(gene["gene_id"])
            if hit is None:
                continue
            retained += 1
            output_handle.write(
                "\t".join(
                    (
                        *(gene[column] for column in _GENE_COLUMNS),
                        hit.subject,
                        hit.identity,
                        hit.alignment_length,
                        hit.query_coverage,
                        hit.taxid,
                    )
                )
                + "\n"
            )

    print(
        f"Kept the best hit of {retained} of {len(genes)} SSU and LSU genes, "
        f"discarding alignments shorter than {minimum_alignment_length} bases."
    )


if "snakemake" in globals():
    workflow = globals()["snakemake"]
    log_path = Path(workflow.log[0])
    with log_path.open("w", encoding="utf-8") as log_handle:
        with redirect_stdout(log_handle), redirect_stderr(log_handle):
            try:
                select_rrna_best_hits(
                    Path(workflow.input.genes),
                    Path(workflow.input.ssu),
                    Path(workflow.input.lsu),
                    Path(workflow.output.hits),
                    minimum_alignment_length=int(workflow.params.minimum),
                )
            except Exception:
                traceback.print_exc()
                raise
