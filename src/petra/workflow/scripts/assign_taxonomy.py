import math
import traceback
from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

TAXONOMIC_RANKS = (
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
_MISSING_TAXON = "N/A"
_REASON_COLUMN = "unassigned_reason"
_NO_REASON = "-"
_INSUFFICIENT_SUPPORT_REASON = "insufficient_support"
_NO_CONSENSUS_REASON = "no_consensus"


@dataclass
class AssemblySequence:
    length: int = 0
    gc_bases: int = 0



@dataclass
class ProteinAssignment:
    bitscore: float
    lineages: list[tuple[str, ...]]


def _fasta_identifiers(source: Path) -> Iterator[str]:
    with source.open(encoding="utf-8") as input_handle:
        for line_number, line in enumerate(input_handle, start=1):
            if not line.startswith(">"):
                continue
            identifier = line[1:].split(maxsplit=1)[0]
            if not identifier:
                raise ValueError(
                    f"Invalid FASTA header at line {line_number} in {source}: "
                    "expected a sequence identifier."
                )
            yield identifier


def _assembly_sequence_id(protein_id: str) -> str:
    sequence_id, separator, gene_number = protein_id.rpartition("_")
    if not separator or not sequence_id or not gene_number.isdecimal():
        raise ValueError(f"Unexpected Pyrodigal protein identifier: {protein_id!r}")
    return sequence_id


def _assembly_sequences(source: Path) -> dict[str, AssemblySequence]:
    sequences: dict[str, AssemblySequence] = {}
    current_sequence: AssemblySequence | None = None
    with source.open(encoding="utf-8") as input_handle:
        for line_number, line in enumerate(input_handle, start=1):
            if line.startswith(">"):
                sequence_id = line[1:].split(maxsplit=1)[0]
                if not sequence_id:
                    raise ValueError(
                        f"Invalid FASTA header at line {line_number} in {source}: "
                        "expected a sequence identifier."
                    )
                if sequence_id in sequences:
                    raise ValueError(
                        f"Duplicate assembly sequence identifier: {sequence_id!r}."
                    )
                current_sequence = AssemblySequence()
                sequences[sequence_id] = current_sequence
            elif line.strip():
                if current_sequence is None:
                    raise ValueError(
                        f"Invalid assembly FASTA at line {line_number}: "
                        "sequence data occurs before the first header."
                    )
                sequence = line.strip()
                current_sequence.length += len(sequence)
                current_sequence.gc_bases += sum(base in "GCgc" for base in sequence)
    return sequences


def _unigene_proteins(source: Path, assembly_ids: set[str]) -> dict[str, str]:
    protein_sequences: dict[str, str] = {}
    for protein_id in _fasta_identifiers(source):
        if protein_id in protein_sequences:
            raise ValueError(f"Duplicate unigene protein identifier: {protein_id!r}.")
        sequence_id = _assembly_sequence_id(protein_id)
        if sequence_id not in assembly_ids:
            raise ValueError(
                f"Unigene protein {protein_id!r} belongs to unknown assembly "
                f"sequence {sequence_id!r}."
            )
        protein_sequences[protein_id] = sequence_id
    return protein_sequences


def _alignment_assignments(
    source: Path, protein_sequences: dict[str, str]
) -> dict[str, ProteinAssignment]:
    assignments: dict[str, ProteinAssignment] = {}
    with source.open(encoding="utf-8") as input_handle:
        for line_number, line in enumerate(input_handle, start=1):
            columns = line.rstrip("\n").split("\t")
            if len(columns) != 5:
                raise ValueError(
                    f"Invalid lineage alignment at line {line_number}: expected "
                    "five tab-separated columns."
                )
            protein_id = columns[0]
            if protein_id not in protein_sequences:
                raise ValueError(
                    f"Lineage alignment at line {line_number} refers to unknown "
                    f"unigene protein {protein_id!r}."
                )
            try:
                bitscore = float(columns[3])
            except ValueError as error:
                raise ValueError(
                    f"Invalid bit score at line {line_number}: {columns[3]!r}."
                ) from error
            if not math.isfinite(bitscore):
                raise ValueError(
                    f"Invalid bit score at line {line_number}: {columns[3]!r}."
                )

            lineage = tuple(columns[4].split(";"))
            if len(lineage) != len(TAXONOMIC_RANKS) or any(
                not taxon for taxon in lineage
            ):
                raise ValueError(
                    f"Invalid lineage at line {line_number}: expected "
                    f"{len(TAXONOMIC_RANKS)} non-empty semicolon-separated ranks."
                )

            assignment = assignments.get(protein_id)
            if assignment is None or bitscore > assignment.bitscore:
                assignments[protein_id] = ProteinAssignment(bitscore, [lineage])
            elif bitscore == assignment.bitscore:
                assignment.lineages.append(lineage)
    return assignments


def _protein_lineage(assignment: ProteinAssignment) -> tuple[str, ...]:
    lineage = []
    for index in range(len(TAXONOMIC_RANKS)):
        values = {candidate[index] for candidate in assignment.lineages}
        lineage.append(values.pop() if len(values) == 1 else _MISSING_TAXON)
    return tuple(lineage)


def _write_gene_taxonomy(
    destination: Path,
    sequence_ids: list[str],
    protein_sequences: dict[str, str],
    protein_lineages: dict[str, tuple[str, ...]],
) -> None:
    assigned_proteins: dict[str, list[tuple[str, tuple[str, ...]]]] = {}
    for protein_id, sequence_id in protein_sequences.items():
        lineage = protein_lineages.get(protein_id)
        if lineage is not None:
            assigned_proteins.setdefault(sequence_id, []).append((protein_id, lineage))

    with destination.open("w", encoding="utf-8") as output_handle:
        for sequence_id in sequence_ids:
            output_handle.write(f"{sequence_id}\n")
            for protein_id, lineage in assigned_proteins.get(sequence_id, []):
                assigned_lineage = ";".join(
                    taxon for taxon in lineage if taxon != _MISSING_TAXON
                )
                output_handle.write(
                    f"\t{protein_id}\t{assigned_lineage or _MISSING_TAXON}\n"
                )


def _minimum_aligned_unigenes(
    length: int,
    minimum: int,
    support_length: int,
    unigenes_per_decade: int,
) -> int:
    """Sequence-level evidence floor, increasing logarithmically above the anchor."""
    if length <= support_length:
        return minimum
    return minimum + math.ceil(
        unigenes_per_decade * math.log10(length / support_length)
    )


def _majority_taxon(
    votes: Counter[str], total: int, alpha: float
) -> str:
    """Return the unique taxon whose observed vote share exceeds alpha."""
    winner = _MISSING_TAXON
    highest_count = 0
    tied = False
    for taxon, count in votes.items():
        if taxon == _MISSING_TAXON:
            continue
        if count > highest_count:
            winner = taxon
            highest_count = count
            tied = False
        elif count == highest_count:
            tied = True
    if tied or highest_count / total <= alpha:
        return _MISSING_TAXON
    return winner


def _consensus_lineage(lineages: list[tuple[str, ...]], alpha: float) -> list[str]:
    """Vote rank by rank, keeping only unigenes that agree with every higher rank.

    Independent per-rank votes can pick a descendant whose ancestor lost a tie,
    or one from another lineage (a virus class under Metazoa). A rank without
    consensus leaves no unigene to vote below it; a rank absent from the winning
    clade (N/A in all its lineages) does not interrupt the descent.
    """
    total = len(lineages)
    consensus: list[str] = []
    for rank_index in range(len(TAXONOMIC_RANKS)):
        taxon = _majority_taxon(
            Counter(lineage[rank_index] for lineage in lineages), total, alpha
        )
        consensus.append(taxon)
        lineages = [lineage for lineage in lineages if lineage[rank_index] == taxon]
    return consensus


def assign_taxonomy(
    assembly: Path,
    proteins: Path,
    alignment: Path,
    destination: Path,
    reasons_destination: Path,
    gene_taxonomy_destination: Path,
    alpha: float = 0.5,
    *,
    min_aligned_unigenes: int = 3,
    unigene_support_length: int = 50_000,
    unigenes_per_decade: int = 3,
) -> None:
    if not 0 < alpha < 1:
        raise ValueError(
            "The alpha-majority threshold must be greater than 0 and less than 1."
        )
    if min_aligned_unigenes < 3:
        raise ValueError("The minimum aligned unigene count must be at least 3.")
    if unigene_support_length < 1:
        raise ValueError("The unigene support length must be at least 1 bp.")
    if unigenes_per_decade < 1:
        raise ValueError("The unigene count increase per decade must be at least 1.")

    sequences = _assembly_sequences(assembly)
    sequence_ids = list(sequences)
    protein_sequences = _unigene_proteins(proteins, set(sequences))
    protein_assignments = _alignment_assignments(alignment, protein_sequences)
    protein_lineages = {
        protein_id: _protein_lineage(assignment)
        for protein_id, assignment in protein_assignments.items()
    }
    _write_gene_taxonomy(
        gene_taxonomy_destination,
        sequence_ids,
        protein_sequences,
        protein_lineages,
    )

    sequence_lineages: dict[str, list[tuple[str, ...]]] = {}
    for protein_id, lineage in protein_lineages.items():
        sequence_lineages.setdefault(protein_sequences[protein_id], []).append(lineage)

    assigned_sequences = 0
    with (
        destination.open("w", encoding="utf-8") as output_handle,
        reasons_destination.open("w", encoding="utf-8") as reasons_handle,
    ):
        output_handle.write(
            "sequence_name\tsequence_length\tgc_percent\t"
            "number_of_aligned_unigenes\t"
            + "\t".join(TAXONOMIC_RANKS)
            + "\n"
        )
        reasons_handle.write(f"sequence_name\t{_REASON_COLUMN}\n")
        for sequence_id in sequence_ids:
            lineages = sequence_lineages.get(sequence_id, [])
            total = len(lineages)
            sequence = sequences[sequence_id]
            reason = _NO_REASON
            lineage = [_MISSING_TAXON] * len(TAXONOMIC_RANKS)
            has_taxonomic_evidence = any(
                taxon != _MISSING_TAXON for protein in lineages for taxon in protein
            )
            if has_taxonomic_evidence:
                if total < _minimum_aligned_unigenes(
                    sequence.length,
                    min_aligned_unigenes,
                    unigene_support_length,
                    unigenes_per_decade,
                ):
                    reason = _INSUFFICIENT_SUPPORT_REASON
                else:
                    lineage = _consensus_lineage(lineages, alpha)
                    if all(taxon == _MISSING_TAXON for taxon in lineage):
                        reason = _NO_CONSENSUS_REASON
            if any(taxon != _MISSING_TAXON for taxon in lineage):
                assigned_sequences += 1
            gc_percent = (
                100 * sequence.gc_bases / sequence.length if sequence.length else 0
            )
            output_handle.write(
                f"{sequence_id}\t{sequence.length}\t{gc_percent:.2f}\t{total}\t"
                + "\t".join(lineage)
                + "\n"
            )
            reasons_handle.write(f"{sequence_id}\t{reason}\n")

    print(
        f"Assigned taxonomy to {assigned_sequences} of {len(sequence_ids)} assembly "
        f"sequences from {len(protein_assignments)} aligned unigenes."
    )


if "snakemake" in globals():
    workflow = globals()["snakemake"]
    log_path = Path(workflow.log[0])
    with log_path.open("w", encoding="utf-8") as log_handle:
        with redirect_stdout(log_handle), redirect_stderr(log_handle):
            try:
                assign_taxonomy(
                    Path(workflow.input.assembly),
                    Path(workflow.input.proteins),
                    Path(workflow.input.alignment),
                    Path(workflow.output.taxonomy),
                    Path(workflow.output.reasons),
                    Path(workflow.output.gene_taxonomy),
                    float(workflow.params.alpha),
                    min_aligned_unigenes=int(workflow.params.min_aligned_unigenes),
                    unigene_support_length=int(workflow.params.unigene_support_length),
                    unigenes_per_decade=int(workflow.params.unigenes_per_decade),
                )
            except Exception:
                traceback.print_exc()
                raise
