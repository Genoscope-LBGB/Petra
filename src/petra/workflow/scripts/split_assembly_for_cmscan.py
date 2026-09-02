import heapq
import traceback
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass
from pathlib import Path

_COPY_BLOCK_SIZE = 1 << 20


@dataclass(frozen=True)
class _Record:
    """Byte range of one FASTA record and the number of bases it holds."""

    start: int
    end: int
    bases: int


def _index_records(assembly: Path) -> list[_Record]:
    """Locate every record without holding the assembly in memory."""
    records: list[_Record] = []
    start: int | None = None
    bases = 0
    offset = 0
    with assembly.open("rb") as input_handle:
        for line_number, line in enumerate(input_handle, start=1):
            if line.startswith(b">"):
                if start is not None:
                    records.append(_Record(start, offset, bases))
                start = offset
                bases = 0
            elif start is None:
                if line.strip():
                    raise ValueError(
                        f"Invalid assembly FASTA at line {line_number} in "
                        f"{assembly}: sequence data occurs before the first header."
                    )
            else:
                bases += len(line.strip())
            offset += len(line)
        if start is not None:
            records.append(_Record(start, offset, bases))
    if not records:
        raise ValueError(f"No FASTA record found in {assembly}.")
    return records


def _assign_records(records: list[_Record], chunk_count: int) -> list[list[_Record]]:
    """Balance chunks by bases, longest sequence first, keeping sequences whole.

    cmscan runtime tracks the bases it scans, and record counts are useless as a
    proxy: one chromosome-sized sequence would decide the wall time of the whole
    scatter. Sequences stay whole because the hit coordinates are reported
    against them.
    """
    chunks: list[list[_Record]] = [[] for _ in range(chunk_count)]
    loads = [(0, index) for index in range(chunk_count)]
    heapq.heapify(loads)
    for record in sorted(records, key=lambda record: record.bases, reverse=True):
        bases, index = heapq.heappop(loads)
        chunks[index].append(record)
        heapq.heappush(loads, (bases + record.bases, index))
    return [sorted(chunk, key=lambda record: record.start) for chunk in chunks]


def _copy_record(source, destination, record: _Record) -> None:
    source.seek(record.start)
    remaining = record.end - record.start
    last_byte = b""
    while remaining:
        block = source.read(min(remaining, _COPY_BLOCK_SIZE))
        if not block:
            raise ValueError(
                f"Assembly ended inside the record at byte offset {record.start}."
            )
        destination.write(block)
        remaining -= len(block)
        last_byte = block[-1:]
    if last_byte != b"\n":
        destination.write(b"\n")


def split_assembly_for_cmscan(assembly: Path, chunks: list[Path]) -> None:
    if not chunks:
        raise ValueError("At least one output chunk is required.")

    records = _index_records(assembly)
    assignments = _assign_records(records, len(chunks))
    with assembly.open("rb") as input_handle:
        for destination, assigned in zip(chunks, assignments, strict=True):
            with destination.open("wb") as output_handle:
                for record in assigned:
                    _copy_record(input_handle, output_handle, record)

    bases_per_chunk = [sum(record.bases for record in chunk) for chunk in assignments]
    print(
        f"Split {len(records)} assembly sequences "
        f"({sum(record.bases for record in records)} bases) into {len(chunks)} "
        f"chunks of {min(bases_per_chunk)} to {max(bases_per_chunk)} bases."
    )


if "snakemake" in globals():
    workflow = globals()["snakemake"]
    log_path = Path(workflow.log[0])
    with log_path.open("w", encoding="utf-8") as log_handle:
        with redirect_stdout(log_handle), redirect_stderr(log_handle):
            try:
                split_assembly_for_cmscan(
                    Path(workflow.input.assembly),
                    [Path(chunk) for chunk in workflow.output.chunks],
                )
            except Exception:
                traceback.print_exc()
                raise
