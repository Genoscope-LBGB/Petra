import csv
import traceback
from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path


def _weighted_median(depth_counts: Counter[int], base_count: int) -> str:
    lower_rank = (base_count - 1) // 2
    upper_rank = base_count // 2
    cumulative = 0
    lower_depth: int | None = None

    for depth in sorted(depth_counts):
        cumulative += depth_counts[depth]
        if lower_depth is None and cumulative > lower_rank:
            lower_depth = depth
        if cumulative > upper_rank:
            assert lower_depth is not None
            total = lower_depth + depth
            return str(total // 2) if total % 2 == 0 else f"{total / 2:g}"

    raise ValueError("Cannot calculate median coverage without covered bases.")


def calculate_median_coverage(source: Path, destination: Path) -> None:
    sequence_count = 0
    current_sequence: str | None = None
    depth_counts: Counter[int] = Counter()
    base_count = 0
    completed_sequences: set[str] = set()

    def write_sequence(writer: csv.writer) -> None:
        nonlocal sequence_count
        if current_sequence is None:
            return
        writer.writerow((current_sequence, _weighted_median(depth_counts, base_count)))
        sequence_count += 1

    with (
        source.open(encoding="utf-8", newline="") as input_handle,
        destination.open("w", encoding="utf-8", newline="") as output_handle,
    ):
        reader = csv.reader(input_handle, delimiter="\t")
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("sequence_name", "median_coverage"))

        for line_number, row in enumerate(reader, start=1):
            if len(row) != 4 or not row[0]:
                raise ValueError(
                    f"Invalid coverage interval at line {line_number} in {source}: "
                    "expected four tab-separated columns."
                )
            sequence, start_text, end_text, depth_text = row
            try:
                start = int(start_text)
                end = int(end_text)
                depth = int(depth_text)
            except ValueError as error:
                raise ValueError(
                    f"Invalid numeric coverage value at line {line_number} in {source}."
                ) from error
            if start < 0 or end <= start or depth < 0:
                raise ValueError(
                    f"Invalid coverage interval at line {line_number} in {source}."
                )

            if sequence != current_sequence:
                write_sequence(writer)
                if sequence in completed_sequences:
                    raise ValueError(
                        f"Coverage intervals for sequence {sequence!r} are not contiguous."
                    )
                if current_sequence is not None:
                    completed_sequences.add(current_sequence)
                current_sequence = sequence
                depth_counts.clear()
                base_count = 0

            interval_length = end - start
            depth_counts[depth] += interval_length
            base_count += interval_length

        write_sequence(writer)

    print(f"Calculated median coverage for {sequence_count} assembly sequences.")


if "snakemake" in globals():
    workflow = globals()["snakemake"]
    log_path = Path(workflow.log[0])
    with log_path.open("w", encoding="utf-8") as log_handle:
        with redirect_stdout(log_handle), redirect_stderr(log_handle):
            try:
                calculate_median_coverage(
                    Path(workflow.input.coverage), Path(workflow.output.coverage)
                )
            except Exception:
                traceback.print_exc()
                raise
