import csv
import traceback
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from assign_taxonomy import TAXONOMIC_RANKS

_MISSING_TAXON = "N/A"
_NO_REASON = "-"
_WHOKARYOTE_REASON = "whokaryote"
_METADATA_COLUMNS = (
    "sequence_name",
    "sequence_length",
    "gc_percent",
    "number_of_aligned_unigenes",
)
_TAXONOMY_COLUMNS = (*_METADATA_COLUMNS, *TAXONOMIC_RANKS)
_COVERAGE_COLUMNS = ("sequence_name", "median_coverage")
_REASON_COLUMNS = ("sequence_name", "unassigned_reason")


def _reasons(source: Path) -> dict[str, str]:
    reasons: dict[str, str] = {}
    with source.open(encoding="utf-8", newline="") as input_handle:
        reader = csv.reader(input_handle, delimiter="\t")
        header = tuple(next(reader, ()))
        if header != _REASON_COLUMNS:
            raise ValueError(
                f"Invalid reason header in {source}: expected {_REASON_COLUMNS!r}, "
                f"found {header!r}."
            )
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(_REASON_COLUMNS) or not all(row):
                raise ValueError(
                    f"Invalid reason at line {line_number} in {source}: expected "
                    "two non-empty tab-separated columns."
                )
            sequence_id, reason = row
            if sequence_id in reasons:
                raise ValueError(
                    f"Duplicate reason for sequence {sequence_id!r} at line "
                    f"{line_number} in {source}."
                )
            reasons[sequence_id] = reason
    return reasons


def _median_coverages(source: Path) -> dict[str, str]:
    coverages: dict[str, str] = {}
    with source.open(encoding="utf-8", newline="") as input_handle:
        reader = csv.reader(input_handle, delimiter="\t")
        header = tuple(next(reader, ()))
        if header != _COVERAGE_COLUMNS:
            raise ValueError(
                f"Invalid coverage header in {source}: expected "
                f"{_COVERAGE_COLUMNS!r}, found {header!r}."
            )
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(_COVERAGE_COLUMNS) or not all(row):
                raise ValueError(
                    f"Invalid median coverage at line {line_number} in {source}: "
                    "expected two non-empty tab-separated columns."
                )
            sequence_id, median_coverage = row
            if sequence_id in coverages:
                raise ValueError(
                    f"Duplicate median coverage for sequence {sequence_id!r} at "
                    f"line {line_number} in {source}."
                )
            coverages[sequence_id] = median_coverage
    return coverages


def aggregate_taxonomy(
    raw_taxonomy: Path,
    coverage: Path,
    assignment_reasons: Path,
    whokaryote_reasons: Path,
    destination: Path,
) -> None:
    remaining_coverages = _median_coverages(coverage)
    remaining_assignment_reasons = _reasons(assignment_reasons)
    remaining_whokaryote_reasons = _reasons(whokaryote_reasons)
    seen_sequence_ids: set[str] = set()

    with (
        raw_taxonomy.open(encoding="utf-8", newline="") as input_handle,
        destination.open("w", encoding="utf-8", newline="") as output_handle,
    ):
        reader = csv.reader(input_handle, delimiter="\t")
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        header = tuple(next(reader, ()))
        if header != _TAXONOMY_COLUMNS:
            raise ValueError(
                f"Invalid raw taxonomy header in {raw_taxonomy}: expected "
                f"{_TAXONOMY_COLUMNS!r}, found {header!r}."
            )
        writer.writerow(
            (
                header[0],
                header[1],
                header[2],
                _COVERAGE_COLUMNS[1],
                *header[3:],
                _REASON_COLUMNS[1],
            )
        )

        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(_TAXONOMY_COLUMNS) or not row[0]:
                raise ValueError(
                    f"Invalid raw taxonomy row at line {line_number} in "
                    f"{raw_taxonomy}: expected {len(_TAXONOMY_COLUMNS)} columns "
                    "and a non-empty sequence identifier."
                )
            sequence_id = row[0]
            if sequence_id in seen_sequence_ids:
                raise ValueError(
                    f"Duplicate raw taxonomy sequence identifier {sequence_id!r} "
                    f"at line {line_number} in {raw_taxonomy}."
                )
            seen_sequence_ids.add(sequence_id)

            try:
                median_coverage = remaining_coverages.pop(sequence_id)
            except KeyError as error:
                raise ValueError(
                    f"Median coverages are missing sequence {sequence_id!r}."
                ) from error
            try:
                assignment_reason = remaining_assignment_reasons.pop(sequence_id)
            except KeyError as error:
                raise ValueError(
                    f"Assignment reasons are missing sequence {sequence_id!r}."
                ) from error
            try:
                whokaryote_reason = remaining_whokaryote_reasons.pop(sequence_id)
            except KeyError as error:
                raise ValueError(
                    f"Whokaryote reasons are missing sequence {sequence_id!r}."
                ) from error

            if assignment_reason != _NO_REASON and whokaryote_reason != _NO_REASON:
                raise ValueError(
                    f"Sequence {sequence_id!r} has conflicting unassigned reasons: "
                    f"{assignment_reason!r} and {whokaryote_reason!r}."
                )
            if whokaryote_reason != _NO_REASON:
                if whokaryote_reason != _WHOKARYOTE_REASON:
                    raise ValueError(
                        f"Unexpected Whokaryote reason {whokaryote_reason!r} for "
                        f"sequence {sequence_id!r}."
                    )
                row[len(_METADATA_COLUMNS) :] = [_MISSING_TAXON] * len(
                    TAXONOMIC_RANKS
                )
                reason = whokaryote_reason
            else:
                reason = assignment_reason
            writer.writerow(
                (row[0], row[1], row[2], median_coverage, *row[3:], reason)
            )

    unknown_coverage_ids = remaining_coverages.keys()
    if unknown_coverage_ids:
        formatted_ids = ", ".join(
            repr(identifier) for identifier in sorted(unknown_coverage_ids)
        )
        raise ValueError(
            "Median coverages contain sequence identifiers absent from the raw "
            f"taxonomy table: {formatted_ids}."
        )
    unknown_assignment_ids = remaining_assignment_reasons.keys()
    if unknown_assignment_ids:
        formatted_ids = ", ".join(
            repr(identifier) for identifier in sorted(unknown_assignment_ids)
        )
        raise ValueError(
            "Assignment reasons contain sequence identifiers absent from the raw "
            f"taxonomy table: {formatted_ids}."
        )
    unknown_whokaryote_ids = remaining_whokaryote_reasons.keys()
    if unknown_whokaryote_ids:
        formatted_ids = ", ".join(
            repr(identifier) for identifier in sorted(unknown_whokaryote_ids)
        )
        raise ValueError(
            "Whokaryote reasons contain sequence identifiers absent from the raw "
            f"taxonomy table: {formatted_ids}."
        )

    print(
        "Aggregated sequence metadata, taxonomy, and unassigned reasons for "
        f"{len(seen_sequence_ids)} sequences."
    )


if "snakemake" in globals():
    workflow = globals()["snakemake"]
    log_path = Path(workflow.log[0])
    with log_path.open("w", encoding="utf-8") as log_handle:
        with redirect_stdout(log_handle), redirect_stderr(log_handle):
            try:
                aggregate_taxonomy(
                    Path(workflow.input.raw_taxonomy),
                    Path(workflow.input.coverage),
                    Path(workflow.input.assignment_reasons),
                    Path(workflow.input.whokaryote_reasons),
                    Path(workflow.output.taxonomy),
                )
            except Exception:
                traceback.print_exc()
                raise
