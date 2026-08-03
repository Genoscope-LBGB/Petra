import csv
import traceback
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from assign_taxonomy import TAXONOMIC_RANKS

_MISSING_TAXON = "N/A"
_REASON_COLUMN = "unassigned_reason"
_NO_REASON = "-"
_WHOKARYOTE_REASON = "whokaryote"
_METADATA_COLUMNS = (
    "sequence_name",
    "sequence_length",
    "number_of_aligned_clusters",
    "number_of_genes_in_aligned_clusters",
)
_RAW_COLUMNS = (*_METADATA_COLUMNS, *TAXONOMIC_RANKS)
_WHOKARYOTE_COLUMNS = ("contig", "predicted")
_PETRA_CLASSES = {
    "Eukaryota": "eukaryote",
    "Bacteria": "prokaryote",
    "Archaea": "prokaryote",
    "Viruses": "viral",
    "Viroids": "viral",
}
_WHOKARYOTE_CLASSES = {"eukaryote", "prokaryote"}


def _whokaryote_predictions(source: Path) -> dict[str, str]:
    predictions: dict[str, str] = {}
    with source.open(encoding="utf-8", newline="") as input_handle:
        reader = csv.reader(input_handle, delimiter="\t")
        header = tuple(next(reader, ()))
        if header != _WHOKARYOTE_COLUMNS:
            raise ValueError(
                f"Invalid Whokaryote header in {source}: expected "
                f"{_WHOKARYOTE_COLUMNS!r}, found {header!r}."
            )

        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(_WHOKARYOTE_COLUMNS) or not all(row):
                raise ValueError(
                    f"Invalid Whokaryote prediction at line {line_number} in "
                    f"{source}: expected two non-empty tab-separated columns."
                )
            sequence_id, prediction = row
            if sequence_id in predictions:
                raise ValueError(
                    f"Duplicate Whokaryote sequence identifier {sequence_id!r} "
                    f"at line {line_number} in {source}."
                )
            if prediction not in _WHOKARYOTE_CLASSES:
                raise ValueError(
                    f"Unexpected Whokaryote prediction {prediction!r} for "
                    f"sequence {sequence_id!r} at line {line_number} in {source}."
                )
            predictions[sequence_id] = prediction
    return predictions


def control_taxonomy(
    raw_taxonomy: Path, whokaryote_predictions: Path, reasons_destination: Path
) -> None:
    predictions = _whokaryote_predictions(whokaryote_predictions)
    raw_sequence_ids: set[str] = set()
    disagreements = 0

    with (
        raw_taxonomy.open(encoding="utf-8", newline="") as input_handle,
        reasons_destination.open("w", encoding="utf-8", newline="") as output_handle,
    ):
        reader = csv.reader(input_handle, delimiter="\t")
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        header = tuple(next(reader, ()))
        if header != _RAW_COLUMNS:
            raise ValueError(
                f"Invalid raw taxonomy header in {raw_taxonomy}: expected "
                f"{_RAW_COLUMNS!r}, found {header!r}."
            )
        writer.writerow(("sequence_name", _REASON_COLUMN))

        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(_RAW_COLUMNS) or not row[0]:
                raise ValueError(
                    f"Invalid raw taxonomy row at line {line_number} in "
                    f"{raw_taxonomy}: expected {len(_RAW_COLUMNS)} columns and a "
                    "non-empty sequence identifier."
                )

            sequence_id = row[0]
            if sequence_id in raw_sequence_ids:
                raise ValueError(
                    f"Duplicate raw taxonomy sequence identifier {sequence_id!r} "
                    f"at line {line_number} in {raw_taxonomy}."
                )
            raw_sequence_ids.add(sequence_id)

            domain = row[len(_METADATA_COLUMNS)]
            if domain == _MISSING_TAXON:
                petra_class = None
            else:
                try:
                    petra_class = _PETRA_CLASSES[domain]
                except KeyError as error:
                    raise ValueError(
                        f"Unexpected PETRA domain {domain!r} for sequence "
                        f"{sequence_id!r} at line {line_number} in {raw_taxonomy}."
                    ) from error

            prediction = predictions.get(sequence_id)
            reason = _NO_REASON
            if (
                petra_class is not None
                and prediction is not None
                and petra_class != prediction
            ):
                reason = _WHOKARYOTE_REASON
                disagreements += 1
            writer.writerow((sequence_id, reason))

    unknown_sequence_ids = predictions.keys() - raw_sequence_ids
    if unknown_sequence_ids:
        formatted_ids = ", ".join(
            repr(identifier) for identifier in sorted(unknown_sequence_ids)
        )
        raise ValueError(
            "Whokaryote predictions contain sequence identifiers absent from the "
            f"raw taxonomy table: {formatted_ids}."
        )

    print(
        f"Reported {disagreements} of {len(raw_sequence_ids)} taxonomy assignments "
        "that disagree with Whokaryote."
    )


if "snakemake" in globals():
    workflow = globals()["snakemake"]
    log_path = Path(workflow.log[0])
    with log_path.open("w", encoding="utf-8") as log_handle:
        with redirect_stdout(log_handle), redirect_stderr(log_handle):
            try:
                control_taxonomy(
                    Path(workflow.input.raw_taxonomy),
                    Path(workflow.input.whokaryote_predictions),
                    Path(workflow.output.reasons),
                )
            except Exception:
                traceback.print_exc()
                raise
