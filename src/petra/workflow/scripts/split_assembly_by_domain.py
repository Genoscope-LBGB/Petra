import csv
import traceback
from collections import Counter
from contextlib import ExitStack, redirect_stderr, redirect_stdout
from pathlib import Path
from typing import TextIO

from assign_taxonomy import TAXONOMIC_RANKS

_METADATA_COLUMNS = (
    "sequence_name",
    "sequence_length",
    "gc_percent",
    "median_coverage",
    "number_of_aligned_unigenes",
)
_TAXONOMY_COLUMNS = (*_METADATA_COLUMNS, *TAXONOMIC_RANKS, "unassigned_reason")
_OUTPUT_BY_DOMAIN = {
    "Eukaryota": "eukaryota",
    "Bacteria": "bacteria",
    "Viruses": "viruses",
    "Viroids": "viruses",
    "Archaea": "archaea",
    "N/A": "unassigned",
}


def _domain_assignments(source: Path) -> dict[str, str]:
    assignments: dict[str, str] = {}
    with source.open(encoding="utf-8", newline="") as input_handle:
        reader = csv.reader(input_handle, delimiter="\t")
        header = tuple(next(reader, ()))
        if header != _TAXONOMY_COLUMNS:
            raise ValueError(
                f"Invalid final taxonomy header in {source}: expected "
                f"{_TAXONOMY_COLUMNS!r}, found {header!r}."
            )

        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(_TAXONOMY_COLUMNS) or not row[0]:
                raise ValueError(
                    f"Invalid final taxonomy row at line {line_number} in {source}: "
                    f"expected {len(_TAXONOMY_COLUMNS)} columns and a non-empty "
                    "sequence identifier."
                )
            sequence_id = row[0]
            if sequence_id in assignments:
                raise ValueError(
                    f"Duplicate final taxonomy sequence identifier {sequence_id!r} "
                    f"at line {line_number} in {source}."
                )
            try:
                assignments[sequence_id] = _OUTPUT_BY_DOMAIN[
                    row[len(_METADATA_COLUMNS)]
                ]
            except KeyError as error:
                raise ValueError(
                    f"Unexpected final domain {row[len(_METADATA_COLUMNS)]!r} for "
                    f"sequence {sequence_id!r} at line {line_number} in {source}."
                ) from error
    return assignments


def _open_outputs(stack: ExitStack, destinations: dict[str, Path]) -> dict[str, TextIO]:
    return {
        name: stack.enter_context(destination.open("w", encoding="utf-8", newline=""))
        for name, destination in destinations.items()
    }


def split_assembly_by_domain(
    assembly: Path,
    taxonomy: Path,
    eukaryota: Path,
    bacteria: Path,
    viruses: Path,
    archaea: Path,
    unassigned: Path,
) -> None:
    assignments = _domain_assignments(taxonomy)
    destinations = {
        "eukaryota": eukaryota,
        "bacteria": bacteria,
        "viruses": viruses,
        "archaea": archaea,
        "unassigned": unassigned,
    }
    assembly_ids: set[str] = set()
    counts: Counter[str] = Counter()

    with ExitStack() as stack:
        output_handles = _open_outputs(stack, destinations)
        input_handle = stack.enter_context(assembly.open(encoding="utf-8", newline=""))
        output_handle: TextIO | None = None
        for line_number, line in enumerate(input_handle, start=1):
            if line.startswith(">"):
                sequence_id = line[1:].split(maxsplit=1)[0]
                if not sequence_id:
                    raise ValueError(
                        f"Invalid FASTA header at line {line_number} in {assembly}: "
                        "expected a sequence identifier."
                    )
                if sequence_id in assembly_ids:
                    raise ValueError(
                        f"Duplicate assembly sequence identifier {sequence_id!r} "
                        f"at line {line_number} in {assembly}."
                    )
                try:
                    output_name = assignments[sequence_id]
                except KeyError as error:
                    raise ValueError(
                        f"Assembly sequence {sequence_id!r} is absent from the "
                        f"final taxonomy table {taxonomy}."
                    ) from error
                assembly_ids.add(sequence_id)
                counts[output_name] += 1
                output_handle = output_handles[output_name]
                output_handle.write(line)
            elif output_handle is not None:
                output_handle.write(line)
            elif line.strip():
                raise ValueError(
                    f"Invalid assembly FASTA at line {line_number} in {assembly}: "
                    "sequence data occurs before the first header."
                )

    unknown_sequence_ids = assignments.keys() - assembly_ids
    if unknown_sequence_ids:
        formatted_ids = ", ".join(
            repr(identifier) for identifier in sorted(unknown_sequence_ids)
        )
        raise ValueError(
            "Final taxonomy table contains sequence identifiers absent from the "
            f"assembly: {formatted_ids}."
        )

    print(
        f"Split {len(assembly_ids)} assembly sequences by final domain: "
        + ", ".join(f"{name}={counts[name]}" for name in destinations)
        + "."
    )


if "snakemake" in globals():
    workflow = globals()["snakemake"]
    log_path = Path(workflow.log[0])
    with log_path.open("w", encoding="utf-8") as log_handle:
        with redirect_stdout(log_handle), redirect_stderr(log_handle):
            try:
                split_assembly_by_domain(
                    Path(workflow.input.assembly),
                    Path(workflow.input.taxonomy),
                    Path(workflow.output.eukaryota),
                    Path(workflow.output.bacteria),
                    Path(workflow.output.viruses),
                    Path(workflow.output.archaea),
                    Path(workflow.output.unassigned),
                )
            except Exception:
                traceback.print_exc()
                raise
