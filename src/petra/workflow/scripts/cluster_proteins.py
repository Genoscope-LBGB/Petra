import traceback
from contextlib import redirect_stderr, redirect_stdout

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


def _cluster_sizes(source: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    cluster_started = False
    representative_id: str | None = None
    member_count = 0

    def finish_cluster(line_number: int) -> None:
        if not cluster_started:
            return
        if representative_id is None:
            raise ValueError(
                f"CD-HIT cluster ending at line {line_number} in {source} has no "
                "representative."
            )
        if representative_id in sizes:
            raise ValueError(
                f"Duplicate CD-HIT cluster representative {representative_id!r} "
                f"in {source}."
            )
        sizes[representative_id] = member_count

    with source.open(encoding="utf-8") as input_handle:
        for line_number, line in enumerate(input_handle, start=1):
            line = line.rstrip("\n")
            if line.startswith(">Cluster "):
                finish_cluster(line_number - 1)
                cluster_started = True
                representative_id = None
                member_count = 0
                continue
            if not cluster_started or not line:
                raise ValueError(
                    f"Invalid CD-HIT cluster membership at line {line_number} in "
                    f"{source}."
                )

            member_count += 1
            if line.rstrip().endswith("*"):
                if representative_id is not None:
                    raise ValueError(
                        f"Multiple representatives in one CD-HIT cluster at line "
                        f"{line_number} in {source}."
                    )
                _, separator, description = line.partition(", >")
                identifier, terminator, _ = description.partition("...")
                identifier = identifier.split(maxsplit=1)[0]
                if not separator or not terminator or not identifier:
                    raise ValueError(
                        f"Invalid CD-HIT representative at line {line_number} in "
                        f"{source}."
                    )
                representative_id = identifier

    finish_cluster(line_number if cluster_started else 0)
    if not sizes:
        raise ValueError(f"CD-HIT produced no cluster memberships in {source}.")
    return sizes


def _record_identifier(record: list[bytes]) -> str:
    identifier = record[0][1:].split(maxsplit=1)[0]
    if not identifier:
        raise ValueError("Clustered protein FASTA contains an empty identifier.")
    return identifier.decode("utf-8")


def _sample_clusters(
    source: Path,
    cluster_sizes: dict[str, int],
    destination: BinaryIO,
    sizes_destination: TextIO,
    maximum: int,
    sequence_id: str,
) -> None:
    seed = int.from_bytes(hashlib.sha256(sequence_id.encode("utf-8")).digest())
    generator = random.Random(seed)
    reservoir: list[tuple[int, list[bytes], str, int]] = []
    seen_ids: set[str] = set()
    total = 0

    for index, record in enumerate(_fasta_records(source)):
        total = index + 1
        representative_id = _record_identifier(record)
        if representative_id in seen_ids:
            raise ValueError(
                f"Duplicate clustered protein identifier: {representative_id!r}."
            )
        seen_ids.add(representative_id)
        try:
            cluster_size = cluster_sizes[representative_id]
        except KeyError as error:
            raise ValueError(
                f"Clustered protein {representative_id!r} has no CD-HIT membership."
            ) from error
        candidate = (index, record, representative_id, cluster_size)
        if index < maximum:
            reservoir.append(candidate)
        else:
            replacement = generator.randrange(index + 1)
            if replacement < maximum:
                reservoir[replacement] = candidate

    if total == 0:
        raise ValueError(
            f"CD-HIT produced no clusters for assembly sequence {sequence_id!r}."
        )
    if seen_ids != cluster_sizes.keys():
        unknown_ids = cluster_sizes.keys() - seen_ids
        formatted_ids = ", ".join(
            repr(identifier) for identifier in sorted(unknown_ids)
        )
        raise ValueError(
            "CD-HIT memberships contain representatives absent from the clustered "
            f"protein FASTA: {formatted_ids}."
        )
    for _, record, representative_id, cluster_size in sorted(reservoir):
        destination.writelines(record)
        sizes_destination.write(f"{representative_id}\t{cluster_size}\n")


def cluster_proteins(
    source: Path,
    destination: Path,
    sizes_destination: Path,
    workers: int,
    max_clusters: int,
    runner: Callable[[ProteinGroup], Path] = _run_cdhit,
    log: TextIO | None = None,
) -> None:
    if workers < 1:
        raise ValueError("CD-HIT worker count must be positive.")
    if max_clusters < 1:
        raise ValueError("The per-sequence cluster limit must be positive.")

    with TemporaryDirectory(prefix="petra-cdhit-") as temporary_directory:
        directory = Path(temporary_directory)
        with (
            destination.open("wb") as output_handle,
            sizes_destination.open("w", encoding="utf-8") as sizes_handle,
            ThreadPoolExecutor(max_workers=workers) as executor,
        ):
            sizes_handle.write("cluster_representative\tnumber_of_genes_in_cluster\n")
            for batch in _group_batches(source, directory, workers):
                clustered_outputs = executor.map(runner, batch)
                for group, clustered in zip(batch, clustered_outputs, strict=True):
                    if log is not None and group.log.is_file():
                        log.write(f"== Assembly sequence {group.sequence_id} ==\n")
                        with group.log.open(encoding="utf-8") as group_log:
                            for line in group_log:
                                log.write(line)
                        log.write("\n")
                        log.flush()
                    _sample_clusters(
                        clustered,
                        _cluster_sizes(Path(f"{clustered}.clstr")),
                        output_handle,
                        sizes_handle,
                        max_clusters,
                        group.sequence_id,
                    )


if "snakemake" in globals():
    workflow = globals()["snakemake"]
    log_path = Path(workflow.log[0])
    with log_path.open("w", encoding="utf-8") as log_handle:
        with redirect_stdout(log_handle), redirect_stderr(log_handle):
            try:
                cluster_proteins(
                    Path(workflow.input.proteins),
                    Path(workflow.output.proteins),
                    Path(workflow.output.cluster_sizes),
                    int(workflow.threads),
                    int(workflow.params.maximum),
                    log=log_handle,
                )
            except Exception:
                traceback.print_exc()
                raise
