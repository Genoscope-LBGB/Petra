import traceback
from contextlib import redirect_stderr, redirect_stdout

from pathlib import Path
from typing import TextIO


def _write_record(
    destination: TextIO,
    header: str,
    sequence_lines: list[str],
    low_complexity_residues: int,
    maximum: int,
) -> None:
    if low_complexity_residues <= maximum:
        destination.write(header)
        destination.writelines(sequence_lines)


def filter_fasta(source: Path, destination: Path, maximum: int) -> None:
    if maximum < 0:
        raise ValueError("The low-complexity threshold cannot be negative.")

    header: str | None = None
    sequence_lines: list[str] = []
    low_complexity_residues = 0

    with (
        source.open(encoding="utf-8") as input_handle,
        destination.open("w", encoding="utf-8") as output_handle,
    ):
        for line in input_handle:
            if line.startswith(">"):
                if header is not None:
                    _write_record(
                        output_handle,
                        header,
                        sequence_lines,
                        low_complexity_residues,
                        maximum,
                    )
                header = line
                sequence_lines = []
                low_complexity_residues = 0
            elif header is None:
                if line.strip():
                    raise ValueError(
                        "Protein FASTA contains sequence data before a header."
                    )
            else:
                sequence_lines.append(line)
                low_complexity_residues += sum(
                    character.islower() for character in line
                )

        if header is None:
            raise ValueError("Protein FASTA contains no records.")
        _write_record(
            output_handle,
            header,
            sequence_lines,
            low_complexity_residues,
            maximum,
        )


log_path = Path(snakemake.log[0])
with log_path.open("w", encoding="utf-8") as log_handle:
    with redirect_stdout(log_handle), redirect_stderr(log_handle):
        try:
            filter_fasta(
                Path(snakemake.input.proteins),
                Path(snakemake.output.proteins),
                int(snakemake.params.maximum),
            )
        except Exception:
            traceback.print_exc()
            raise
