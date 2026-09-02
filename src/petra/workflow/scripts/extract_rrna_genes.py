import traceback
from collections.abc import Iterator
from contextlib import ExitStack, redirect_stderr, redirect_stdout
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

# Field positions of a cmscan '--fmt 2' tblout row, Infernal user guide.
_MODEL_NAME = 1
_MODEL_ACCESSION = 2
_SEQUENCE_NAME = 3
_CLAN = 5
_SEQUENCE_FROM = 9
_SEQUENCE_TO = 10
_STRAND = 11
_TRUNCATED = 12
_SCORE = 16
_EVALUE = 17
_INCLUSION = 18
_OVERLAP = 19
_MINIMUM_FIELDS = 27

# Hits marked '=' overlap a higher-scoring hit of the same clan. cmscan runs with
# --clanin, so this is where the several per-domain SSU and LSU models covering
# one rRNA region are reduced to their best hit.
_OVERLAPPED = "="

# Rfam families of RRNA_CM_ACCESSIONS that have a BLAST database. 5S (RF00001)
# and 5.8S (RF00002) have none and reach the GFF inventory only.
_SUBUNIT_BY_ACCESSION = {
    "RF00177": "SSU",
    "RF01959": "SSU",
    "RF01960": "SSU",
    "RF02542": "SSU",
    "RF02545": "SSU",
    "RF02540": "LSU",
    "RF02541": "LSU",
    "RF02543": "LSU",
    "RF02546": "LSU",
}
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
_COMPLEMENTS = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "U": "A",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "V": "B",
    "D": "H",
    "H": "D",
    "N": "N",
}
_COMPLEMENT_TABLE = str.maketrans(
    "".join(_COMPLEMENTS) + "".join(_COMPLEMENTS).lower(),
    "".join(_COMPLEMENTS.values()) + "".join(_COMPLEMENTS.values()).lower(),
)


@dataclass(frozen=True)
class RnaGene:
    """One cmscan hit, with coordinates in GFF orientation (start <= end)."""

    gene_id: str
    sequence_name: str
    start: int
    end: int
    strand: str
    accession: str
    family: str
    clan: str
    score: str
    evalue: str
    truncated: str
    inclusion: str
    subunit: str | None


def _subunit(accession: str, family: str) -> str | None:
    subunit = _SUBUNIT_BY_ACCESSION.get(accession)
    if subunit is not None:
        return subunit
    for candidate in ("SSU", "LSU"):
        if family.startswith(f"{candidate}_"):
            return candidate
    return None


def parse_tblout(source: Path) -> list[RnaGene]:
    genes: list[RnaGene] = []
    with source.open(encoding="utf-8") as input_handle:
        for line_number, line in enumerate(input_handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if len(fields) < _MINIMUM_FIELDS:
                raise ValueError(
                    f"Expected at least {_MINIMUM_FIELDS} fields at line "
                    f"{line_number} of {source}, got {len(fields)}. cmscan must "
                    "run with '--fmt 2'."
                )
            strand = fields[_STRAND]
            if strand not in ("+", "-"):
                raise ValueError(
                    f"Unexpected strand {strand!r} at line {line_number} of "
                    f"{source}. cmscan must run with '--fmt 2'."
                )
            if fields[_OVERLAP] == _OVERLAPPED:
                continue
            try:
                first = int(fields[_SEQUENCE_FROM])
                second = int(fields[_SEQUENCE_TO])
            except ValueError as error:
                raise ValueError(
                    f"Invalid hit coordinates at line {line_number} of {source}."
                ) from error
            if first < 1 or second < 1:
                raise ValueError(
                    f"Non-positive hit coordinate at line {line_number} of {source}."
                )
            accession = fields[_MODEL_ACCESSION]
            family = fields[_MODEL_NAME]
            genes.append(
                RnaGene(
                    gene_id=f"rrna_{len(genes) + 1}",
                    sequence_name=fields[_SEQUENCE_NAME],
                    start=min(first, second),
                    end=max(first, second),
                    strand=strand,
                    accession=accession,
                    family=family,
                    clan=fields[_CLAN],
                    score=fields[_SCORE],
                    evalue=fields[_EVALUE],
                    truncated=fields[_TRUNCATED],
                    inclusion=fields[_INCLUSION],
                    subunit=_subunit(accession, family),
                )
            )
    return genes


def _write_gff(genes: list[RnaGene], destination: Path) -> None:
    with destination.open("w", encoding="utf-8", newline="\n") as output_handle:
        output_handle.write("##gff-version 3\n")
        for gene in genes:
            attributes = ";".join(
                (
                    f"ID={gene.gene_id}",
                    f"rfam_accession={gene.accession}",
                    f"rfam_family={gene.family}",
                    f"clan={gene.clan}",
                    f"evalue={gene.evalue}",
                    f"truncated={gene.truncated}",
                    f"inclusion={gene.inclusion}",
                )
            )
            output_handle.write(
                "\t".join(
                    (
                        gene.sequence_name,
                        "cmscan",
                        "rRNA",
                        str(gene.start),
                        str(gene.end),
                        gene.score,
                        gene.strand,
                        ".",
                        attributes,
                    )
                )
                + "\n"
            )


def _write_genes(genes: list[RnaGene], destination: Path) -> None:
    with destination.open("w", encoding="utf-8", newline="\n") as output_handle:
        output_handle.write("\t".join(_GENE_COLUMNS) + "\n")
        for gene in genes:
            if gene.subunit is None:
                continue
            output_handle.write(
                "\t".join(
                    (
                        gene.gene_id,
                        gene.sequence_name,
                        str(gene.start),
                        str(gene.end),
                        gene.strand,
                        gene.accession,
                        gene.subunit,
                        gene.score,
                    )
                )
                + "\n"
            )


def _iter_fasta(source: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    chunks: list[str] = []
    with source.open(encoding="utf-8") as input_handle:
        for line_number, line in enumerate(input_handle, start=1):
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].split(maxsplit=1)[0] if line[1:].strip() else ""
                if not header:
                    raise ValueError(
                        f"Invalid FASTA header at line {line_number} in {source}: "
                        "expected a sequence identifier."
                    )
                chunks = []
            elif header is None:
                if line.strip():
                    raise ValueError(
                        f"Invalid FASTA at line {line_number} in {source}: "
                        "sequence data occurs before the first header."
                    )
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def _write_gene_sequences(
    assembly: Path,
    genes: list[RnaGene],
    handles: dict[str, TextIO],
) -> int:
    genes_by_sequence: dict[str, list[RnaGene]] = {}
    for gene in genes:
        if gene.subunit is not None:
            genes_by_sequence.setdefault(gene.sequence_name, []).append(gene)

    written = 0
    for sequence_name, sequence in _iter_fasta(assembly):
        pending = genes_by_sequence.pop(sequence_name, None)
        if pending is None:
            continue
        for gene in pending:
            if gene.end > len(sequence):
                raise ValueError(
                    f"Hit {gene.gene_id} ends at {gene.end} but sequence "
                    f"{sequence_name} is {len(sequence)} bases long."
                )
            subsequence = sequence[gene.start - 1 : gene.end]
            if gene.strand == "-":
                subsequence = subsequence.translate(_COMPLEMENT_TABLE)[::-1]
            handles[gene.subunit].write(f">{gene.gene_id}\n{subsequence}\n")
            written += 1

    if genes_by_sequence:
        missing = ", ".join(repr(name) for name in sorted(genes_by_sequence))
        raise ValueError(
            f"cmscan reported hits on sequences absent from {assembly}: {missing}."
        )
    return written


def extract_rrna_genes(
    assembly: Path,
    tblout: Path,
    gff: Path,
    genes_table: Path,
    ssu: Path,
    lsu: Path,
) -> None:
    genes = parse_tblout(tblout)
    _write_gff(genes, gff)
    _write_genes(genes, genes_table)

    with ExitStack() as stack:
        handles = {
            "SSU": stack.enter_context(ssu.open("w", encoding="utf-8", newline="\n")),
            "LSU": stack.enter_context(lsu.open("w", encoding="utf-8", newline="\n")),
        }
        written = _write_gene_sequences(assembly, genes, handles)

    unclassifiable = sum(1 for gene in genes if gene.subunit is None)
    if not genes:
        print(f"No rRNA gene detected in {assembly}.")
    elif not written:
        print(
            f"Detected {len(genes)} rRNA genes, none of them SSU or LSU: "
            "nothing to classify."
        )
    else:
        print(
            f"Detected {len(genes)} rRNA genes, extracted {written} SSU and LSU "
            f"genes, left {unclassifiable} genes without a BLAST database."
        )


if "snakemake" in globals():
    workflow = globals()["snakemake"]
    log_path = Path(workflow.log[0])
    with log_path.open("w", encoding="utf-8") as log_handle:
        with redirect_stdout(log_handle), redirect_stderr(log_handle):
            try:
                extract_rrna_genes(
                    Path(workflow.input.assembly),
                    Path(workflow.input.tblout),
                    Path(workflow.output.gff),
                    Path(workflow.output.genes),
                    Path(workflow.output.ssu),
                    Path(workflow.output.lsu),
                )
            except Exception:
                traceback.print_exc()
                raise
