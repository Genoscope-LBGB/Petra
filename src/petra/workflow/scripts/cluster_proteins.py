import hashlib
import random
import subprocess
from collections.abc import Callable, Iterator
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import BinaryIO, TextIO


@dataclass(frozen=True)
class ProteinGroup:
    sequence_id: str
    source: Path
    clustered: Path
    log: Path


def _sequence_id(header: str) -> str:
    identifier = header[1:].split(maxsplit=1)[0]
    sequence_id, separator, gene_number = identifier.rpartition("_")
    if not separator or not sequence_id or not gene_number.isdecimal():
        raise ValueError(f"Unexpected Pyrodigal protein identifier: {identifier!r}")
    return sequence_id


def _group_batches(
    source: Path,
    directory: Path,
    batch_size: int,
) -> Iterator[list[ProteinGroup]]:
    batch: list[ProteinGroup] = []
    seen_sequence_ids: set[str] = set()
    current_sequence_id: str | None = None
    current_handle: TextIO | None = None
    group_number = 0

    with source.open(encoding="utf-8") as input_handle:
        try:
            for line in input_handle:
                if line.startswith(">"):
                    sequence_id = _sequence_id(line)
                    if sequence_id != current_sequence_id:
                        if current_handle is not None:
                            current_handle.close()
                            current_handle = None
                        if current_sequence_id is not None:
                            seen_sequence_ids.add(current_sequence_id)
                        if len(batch) == batch_size:
                            yield batch
                            batch = []
                        if sequence_id in seen_sequence_ids:
                            raise ValueError(
                                "Proteins from assembly sequence "
                                f"{sequence_id!r} are not contiguous."
                            )

                        group_number += 1
                        prefix = directory / f"{group_number:08d}"
                        group = ProteinGroup(
                            sequence_id=sequence_id,
                            source=prefix.with_suffix(".faa"),
                            clustered=prefix.with_suffix(".clustered.faa"),
                            log=prefix.with_suffix(".log"),
                        )
                        batch.append(group)
                        current_handle = group.source.open("w", encoding="utf-8")
                        current_sequence_id = sequence_id
                    current_handle.write(line)
                elif current_handle is None:
                    if line.strip():
                        raise ValueError(
                            "Protein FASTA contains sequence data before a header."
                        )
                else:
                    current_handle.write(line)
        finally:
            if current_handle is not None:
                current_handle.close()

    if batch:
        yield batch


def _run_cdhit(group: ProteinGroup) -> Path:
    command = [
        "cd-hit",
        "-i",
        str(group.source),
        "-o",
        str(group.clustered),
        "-n",
        "5",
        "-d",
        "0",
        "-T",
        "1",
    ]
    with group.log.open("w", encoding="utf-8") as log_handle:
        result = subprocess.run(
            command,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if result.returncode != 0:
        details = group.log.read_text(encoding="utf-8", errors="replace")
        raise RuntimeError(
            f"CD-HIT failed for assembly sequence {group.sequence_id!r}:\n{details}"
        )
    return group.clustered


def _fasta_records(source: Path) -> Iterator[list[bytes]]:
    record: list[bytes] | None = None
    with source.open("rb") as input_handle:
        for line in input_handle:
            if line.startswith(b">"):
                if record is not None:
                    yield record
                record = [line]
            elif record is None:
                if line.strip():
                    raise ValueError(
                        "Clustered protein FASTA contains sequence data before a header."
                    )
            else:
                record.append(line)
    if record is not None:
        yield record


def _sample_clusters(
    source: Path,
    destination: BinaryIO,
    maximum: int,
    sequence_id: str,
) -> None:
    seed = int.from_bytes(hashlib.sha256(sequence_id.encode("utf-8")).digest())
    generator = random.Random(seed)
    reservoir: list[tuple[int, list[bytes]]] = []
    total = 0

    for index, record in enumerate(_fasta_records(source)):
        total = index + 1
        if index < maximum:
            reservoir.append((index, record))
        else:
            replacement = generator.randrange(index + 1)
            if replacement < maximum:
                reservoir[replacement] = (index, record)

    if total == 0:
        raise ValueError(
            f"CD-HIT produced no clusters for assembly sequence {sequence_id!r}."
        )
    for _, record in sorted(reservoir):
        destination.writelines(record)


def cluster_proteins(
    source: Path,
    destination: Path,
    workers: int,
    max_clusters: int,
    runner: Callable[[ProteinGroup], Path] = _run_cdhit,
) -> None:
    if workers < 1:
        raise ValueError("CD-HIT worker count must be positive.")
    if max_clusters < 1:
        raise ValueError("The per-sequence cluster limit must be positive.")

    with TemporaryDirectory(prefix="petra-cdhit-") as temporary_directory:
        directory = Path(temporary_directory)
        with destination.open("wb") as output_handle, ThreadPoolExecutor(
            max_workers=workers
        ) as executor:
            for batch in _group_batches(source, directory, workers):
                clustered_outputs = executor.map(runner, batch)
                for group, clustered in zip(
                    batch, clustered_outputs, strict=True
                ):
                    _sample_clusters(
                        clustered,
                        output_handle,
                        max_clusters,
                        group.sequence_id,
                    )


if "snakemake" in globals():
    cluster_proteins(
        Path(snakemake.input.proteins),
        Path(snakemake.output.proteins),
        int(snakemake.threads),
        int(snakemake.params.maximum),
    )
