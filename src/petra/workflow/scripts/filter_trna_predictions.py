import traceback
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

_COLUMNS = (
    "sequence_name",
    "trna_number",
    "start",
    "end",
    "type",
    "codon",
    "intron_start",
    "intron_end",
    "score",
    "note",
)
# tRNAscan-SE prints three header lines before its tabular predictions.
_HEADER_LINES = 3
_HEADER_PREFIXES = ("Sequence", "Name", "-")
_MINIMUM_FIELDS = 9
_TYPE = 4
_NOTE = 9
_NO_NOTE = "-"
_PSEUDOGENE = "pseudo"
_EXCLUDED_TYPES = frozenset({"Undet", "Sup"})


def filter_trna_predictions(source: Path, destination: Path) -> None:
    """Keep the confident tRNA predictions and give them named columns.

    Pseudogenes and the undetermined and suppressor isotypes are the predictions
    tRNAscan-SE itself marks as unreliable or non-standard.
    """
    kept = 0
    discarded = 0
    with (
        source.open(encoding="utf-8") as input_handle,
        destination.open("w", encoding="utf-8", newline="\n") as output_handle,
    ):
        output_handle.write("\t".join(_COLUMNS) + "\n")
        for line_number, line in enumerate(input_handle, start=1):
            if line_number <= _HEADER_LINES:
                if not line.startswith(_HEADER_PREFIXES):
                    raise ValueError(
                        f"Unexpected tRNAscan-SE header at line {line_number} of "
                        f"{source}: {line.strip()!r}."
                    )
                continue
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) < _MINIMUM_FIELDS:
                raise ValueError(
                    f"Expected at least {_MINIMUM_FIELDS} fields at line "
                    f"{line_number} of {source}, got {len(fields)}."
                )
            note = fields[_NOTE] if len(fields) > _NOTE else _NO_NOTE
            if note == _PSEUDOGENE or fields[_TYPE] in _EXCLUDED_TYPES:
                discarded += 1
                continue
            kept += 1
            output_handle.write("\t".join((*fields[:_MINIMUM_FIELDS], note)) + "\n")

    print(
        f"Kept {kept} tRNA predictions and discarded {discarded} pseudogene, "
        "undetermined or suppressor predictions."
    )


if "snakemake" in globals():
    workflow = globals()["snakemake"]
    log_path = Path(workflow.log[0])
    with log_path.open("w", encoding="utf-8") as log_handle:
        with redirect_stdout(log_handle), redirect_stderr(log_handle):
            try:
                filter_trna_predictions(
                    Path(workflow.input.predictions),
                    Path(workflow.output.predictions),
                )
            except Exception:
                traceback.print_exc()
                raise
