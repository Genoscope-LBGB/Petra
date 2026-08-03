import csv
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
_GENE_DENSITY_REASON = "gene_density"
_NO_CONSENSUS_REASON = "no_consensus"
_CLUSTER_SIZE_COLUMNS = (
    "cluster_representative",
    "number_of_genes_in_cluster",
)


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


def _assembly_sequences(source: Path) -> dict[str, int]:
    sequence_lengths: dict[str, int] = {}
    current_sequence_id: str | None = None
    with source.open(encoding="utf-8") as input_handle:
        for line_number, line in enumerate(input_handle, start=1):
            if line.startswith(">"):
                sequence_id = line[1:].split(maxsplit=1)[0]
                if not sequence_id:
                    raise ValueError(
                        f"Invalid FASTA header at line {line_number} in {source}: "
                        "expected a sequence identifier."
                    )
                if sequence_id in sequence_lengths:
                    raise ValueError(
                        f"Duplicate assembly sequence identifier: {sequence_id!r}."
                    )
                sequence_lengths[sequence_id] = 0
                current_sequence_id = sequence_id
            elif line.strip():
                if current_sequence_id is None:
                    raise ValueError(
                        f"Invalid assembly FASTA at line {line_number}: "
                        "sequence data occurs before the first header."
                    )
                sequence_lengths[current_sequence_id] += len(line.strip())
    return sequence_lengths


def _clustered_proteins(source: Path, assembly_ids: set[str]) -> dict[str, str]:
    protein_sequences: dict[str, str] = {}
    for protein_id in _fasta_identifiers(source):
        if protein_id in protein_sequences:
            raise ValueError(f"Duplicate clustered protein identifier: {protein_id!r}.")
        sequence_id = _assembly_sequence_id(protein_id)
        if sequence_id not in assembly_ids:
            raise ValueError(
                f"Clustered protein {protein_id!r} belongs to unknown assembly "
                f"sequence {sequence_id!r}."
            )
        protein_sequences[protein_id] = sequence_id
    return protein_sequences


def _cluster_sizes(source: Path, protein_sequences: dict[str, str]) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with source.open(encoding="utf-8", newline="") as input_handle:
        reader = csv.reader(input_handle, delimiter="\t")
        header = tuple(next(reader, ()))
        if header != _CLUSTER_SIZE_COLUMNS:
            raise ValueError(
                f"Invalid cluster-size header in {source}: expected "
                f"{_CLUSTER_SIZE_COLUMNS!r}, found {header!r}."
            )
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(_CLUSTER_SIZE_COLUMNS) or not all(row):
                raise ValueError(
                    f"Invalid cluster size at line {line_number} in {source}: "
                    "expected two non-empty tab-separated columns."
                )
            representative_id, size_text = row
            if representative_id not in protein_sequences:
                raise ValueError(
                    f"Cluster size at line {line_number} refers to unknown clustered "
                    f"protein {representative_id!r}."
                )
            if representative_id in sizes:
                raise ValueError(
                    f"Duplicate cluster size for {representative_id!r} at line "
                    f"{line_number} in {source}."
                )
            try:
                size = int(size_text)
            except ValueError as error:
                raise ValueError(
                    f"Invalid cluster size at line {line_number} in {source}: "
                    f"{size_text!r}."
                ) from error
            if size < 1:
                raise ValueError(
                    f"Invalid cluster size at line {line_number} in {source}: "
                    "expected a positive integer."
                )
            sizes[representative_id] = size

    missing_ids = protein_sequences.keys() - sizes.keys()
    if missing_ids:
        formatted_ids = ", ".join(
            repr(identifier) for identifier in sorted(missing_ids)
        )
        raise ValueError(
            "Clustered proteins are missing cluster sizes: " f"{formatted_ids}."
        )
    return sizes


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
                    f"clustered protein {protein_id!r}."
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


def _majority_taxon(votes: Counter[str], total: int, alpha: float) -> str:
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
    if tied or highest_count / total < alpha:
        return _MISSING_TAXON
    return winner


def assign_taxonomy(
    assembly: Path,
    proteins: Path,
    cluster_sizes: Path,
    alignment: Path,
    destination: Path,
    reasons_destination: Path,
    alpha: float = 0.5,
    max_bases_per_aligned_gene: int = 50_000,
) -> None:
    if not 0 < alpha <= 1:
        raise ValueError(
            "The alpha-majority threshold must be greater than 0 and at most 1."
        )
    if max_bases_per_aligned_gene < 1:
        raise ValueError("The maximum bases per aligned gene must be positive.")

    sequence_lengths = _assembly_sequences(assembly)
    sequence_ids = list(sequence_lengths)
    protein_sequences = _clustered_proteins(proteins, set(sequence_lengths))
    protein_cluster_sizes = _cluster_sizes(cluster_sizes, protein_sequences)
    protein_assignments = _alignment_assignments(alignment, protein_sequences)

    rank_votes: dict[str, list[Counter[str]]] = {}
    assignment_counts: Counter[str] = Counter()
    aligned_gene_counts: Counter[str] = Counter()
    for protein_id, assignment in protein_assignments.items():
        sequence_id = protein_sequences[protein_id]
        assignment_counts[sequence_id] += 1
        aligned_gene_counts[sequence_id] += protein_cluster_sizes[protein_id]
        sequence_votes = rank_votes.setdefault(
            sequence_id, [Counter() for _ in TAXONOMIC_RANKS]
        )
        for votes, taxon in zip(
            sequence_votes, _protein_lineage(assignment), strict=True
        ):
            votes[taxon] += 1

    assigned_sequences = 0
    with (
        destination.open("w", encoding="utf-8") as output_handle,
        reasons_destination.open("w", encoding="utf-8") as reasons_handle,
    ):
        output_handle.write(
            "sequence_name\tnumber_of_aligned_clusters\t"
            "number_of_genes_in_aligned_clusters\t"
            + "\t".join(TAXONOMIC_RANKS)
            + "\n"
        )
        reasons_handle.write(f"sequence_name\t{_REASON_COLUMN}\n")
        for sequence_id in sequence_ids:
            total = assignment_counts[sequence_id]
            aligned_gene_count = aligned_gene_counts[sequence_id]
            meets_density = (
                aligned_gene_count > 0
                and sequence_lengths[sequence_id]
                <= aligned_gene_count * max_bases_per_aligned_gene
            )
            reason = _NO_REASON
            if total > 0 and not meets_density:
                lineage = [_MISSING_TAXON] * len(TAXONOMIC_RANKS)
                reason = _GENE_DENSITY_REASON
            elif meets_density:
                votes_by_rank = rank_votes[sequence_id]
                lineage = [
                    _majority_taxon(votes, total, alpha) for votes in votes_by_rank
                ]
                has_taxonomic_evidence = any(
                    taxon != _MISSING_TAXON
                    for votes in votes_by_rank
                    for taxon in votes
                )
                if has_taxonomic_evidence and all(
                    taxon == _MISSING_TAXON for taxon in lineage
                ):
                    reason = _NO_CONSENSUS_REASON
            else:
                lineage = [_MISSING_TAXON] * len(TAXONOMIC_RANKS)
            if any(taxon != _MISSING_TAXON for taxon in lineage):
                assigned_sequences += 1
            output_handle.write(
                f"{sequence_id}\t{total}\t{aligned_gene_count}\t"
                + "\t".join(lineage)
                + "\n"
            )
            reasons_handle.write(f"{sequence_id}\t{reason}\n")

    print(
        f"Assigned taxonomy to {assigned_sequences} of {len(sequence_ids)} assembly "
        f"sequences from {len(protein_assignments)} aligned clusters representing "
        f"{sum(aligned_gene_counts.values())} genes."
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
                    Path(workflow.input.cluster_sizes),
                    Path(workflow.input.alignment),
                    Path(workflow.output.taxonomy),
                    Path(workflow.output.reasons),
                    float(workflow.params.alpha),
                    int(workflow.params.max_bases_per_aligned_gene),
                )
            except Exception:
                traceback.print_exc()
                raise
