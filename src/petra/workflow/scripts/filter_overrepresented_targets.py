import traceback
from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path


def _target_id(line: str, line_number: int) -> str:
    columns = line.rstrip("\n").split("\t")
    if len(columns) < 2 or not columns[1]:
        raise ValueError(
            f"Invalid DIAMOND alignment at line {line_number}: "
            "expected a non-empty target sequence identifier in column 2."
        )
    return columns[1]


def filter_alignment(source: Path, destination: Path, maximum: int) -> None:
    if maximum < 0:
        raise ValueError("The target occurrence threshold cannot be negative.")

    target_counts: Counter[str] = Counter()
    with source.open(encoding="utf-8") as input_handle:
        for line_number, line in enumerate(input_handle, start=1):
            target_counts[_target_id(line, line_number)] += 1

    retained = 0
    with (
        source.open(encoding="utf-8") as input_handle,
        destination.open("w", encoding="utf-8") as output_handle,
    ):
        for line_number, line in enumerate(input_handle, start=1):
            if target_counts[_target_id(line, line_number)] <= maximum:
                output_handle.write(line)
                retained += 1

    total = sum(target_counts.values())
    removed_targets = sum(count > maximum for count in target_counts.values())
    print(
        f"Retained {retained} of {total} alignments; "
        f"removed {total - retained} alignments to "
        f"{removed_targets} overrepresented target sequences."
    )


if "snakemake" in globals():
    workflow = globals()["snakemake"]
    log_path = Path(workflow.log[0])
    with log_path.open("w", encoding="utf-8") as log_handle:
        with redirect_stdout(log_handle), redirect_stderr(log_handle):
            try:
                filter_alignment(
                    Path(workflow.input.alignment),
                    Path(workflow.output.alignment),
                    int(workflow.params.maximum),
                )
            except Exception:
                traceback.print_exc()
                raise
