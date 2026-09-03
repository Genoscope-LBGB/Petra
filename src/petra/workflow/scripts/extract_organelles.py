import csv
import traceback
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

_INPUT_COLUMNS = ("sequence_id", "class_fst_stage", "class_snd_stage")
_OUTPUT_COLUMNS = ("sequence_id", "organelle_type")
_ORGANELLE_CLASS = "organelle"


def extract_organelles(source: Path, destination: Path) -> None:
    """Write the identifiers and types of Tiara organelle predictions."""
    count = 0
    with (
        source.open(encoding="utf-8", newline="") as input_handle,
        destination.open("w", encoding="utf-8", newline="") as output_handle,
    ):
        reader = csv.reader(input_handle, delimiter="\t")
        header = tuple(next(reader, ()))
        if header != _INPUT_COLUMNS:
            raise ValueError(
                f"Invalid Tiara header in {source}: expected {_INPUT_COLUMNS!r}, "
                f"found {header!r}."
            )

        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(_OUTPUT_COLUMNS)
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(_INPUT_COLUMNS) or not all(row):
                raise ValueError(
                    f"Invalid Tiara prediction at line {line_number} in {source}: "
                    "expected three non-empty tab-separated columns."
                )
            sequence_id, first_stage_class, second_stage_class = row
            if first_stage_class != _ORGANELLE_CLASS:
                continue
            writer.writerow((sequence_id, second_stage_class))
            count += 1

    print(f"Reported {count} Tiara organelle predictions.")


if "snakemake" in globals():
    workflow = globals()["snakemake"]
    log_path = Path(workflow.log[0])
    with (
        log_path.open("w", encoding="utf-8") as log_handle,
        redirect_stdout(log_handle),
        redirect_stderr(log_handle),
    ):
        try:
            extract_organelles(
                Path(workflow.input.predictions),
                Path(workflow.output.organelles),
            )
        except Exception:
            traceback.print_exc()
            raise
