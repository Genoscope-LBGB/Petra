from __future__ import annotations

import logging
import shutil
from importlib.metadata import PackageNotFoundError, version
from enum import Enum
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.logging import RichHandler
from rich.table import Table

from .database import (
    CLAN_INFO_FILENAME,
    COVARIANCE_MODEL_FILENAME,
    DIAMOND_DATABASE_FILENAME,
    LSU_DATABASE_FILENAME,
    SSU_DATABASE_FILENAME,
    DatabaseBuildError,
    build_diamond_database,
    build_rna_databases,
)
from .launcher import LauncherError, launch_workflow, profile_summaries
from .taxonomy import HostTaxonError, resolve_host_taxon

_HELP_CONTEXT = {"help_option_names": ["-h", "--help"]}
_RUN_CONTEXT = {**_HELP_CONTEXT, "allow_extra_args": True}
_RRNA_CHUNKS = 16


class SequencingTechnology(str, Enum):
    PACBIO = "pacbio"
    ONT = "ont"


class SummaryPanel(str, Enum):
    HEADER = "header"
    GC_HISTOGRAM = "gc-histogram"
    COVERAGE_HISTOGRAM = "coverage-histogram"
    TAXONOMY_LEGEND = "taxonomy-legend"
    LONGEST_SEQUENCES = "longest-sequences"
    LENGTH_LEGEND = "length-legend"
    ASSIGNMENT_STATUS = "assignment-status"


def _require_executables(*, rrna: bool, trna: bool) -> None:
    required = [
        "diamond",
        "pyrodigal",
        "taxonkit",
        "cd-hit",
        "minimap2",
        "samtools",
        "segmasker",
        "whokaryote.py",
        "bedtools",
    ]
    if rrna:
        required.extend(("cmscan", "blastn"))
    if trna:
        required.append("tRNAscan-SE")
    missing = [
        executable for executable in required if shutil.which(executable) is None
    ]
    if missing:
        raise LauncherError(
            f"Required executables not found on PATH: {', '.join(missing)}."
        )


def _validate_alpha(value: float) -> float:
    if value <= 0:
        raise typer.BadParameter(
            "The alpha-majority threshold must be greater than 0."
        )
    return value


def _validate_database_directory(value: Path) -> Path:
    if not (value / DIAMOND_DATABASE_FILENAME).is_file():
        raise typer.BadParameter(
            f"{DIAMOND_DATABASE_FILENAME} not found in {value}. Build it with "
            f"'petra build-db --output {value}'."
        )
    return value


app = typer.Typer(
    help="PETRA - Predicted protein-enabled Taxonomic Rank Assignment",
    context_settings=_HELP_CONTEXT,
    no_args_is_help=True,
    pretty_exceptions_show_locals=False,
)
console = Console()


def _configure_logging() -> None:
    petra_logger = logging.getLogger("petra")
    if not petra_logger.handlers:
        handler = RichHandler(
            console=console,
            show_time=False,
            show_path=False,
            markup=False,
        )
        handler.setFormatter(logging.Formatter("%(message)s"))
        petra_logger.addHandler(handler)
    petra_logger.setLevel(logging.INFO)
    petra_logger.propagate = False


def _version_callback(value: bool) -> None:
    if value:
        try:
            workflow_version = version("petra")
        except PackageNotFoundError:
            workflow_version = "0+unknown"
        typer.echo(workflow_version)
        raise typer.Exit()


@app.callback()
def main(
    version_: Annotated[
        bool | None,
        typer.Option(
            "--version",
            callback=_version_callback,
            help="Show the installed workflow version and exit.",
            is_eager=True,
        ),
    ] = None,
) -> None:
    """PETRA - Predicted protein-enabled Taxonomic Rank Assignment"""
    _configure_logging()


@app.command(
    "run",
    context_settings=_RUN_CONTEXT,
    epilog=(
        "Targets and advanced Snakemake options follow '--'. Example: "
        "petra run --assembly assembly.fasta --database database "
        "--output run-1 -- --targets pyrodigal"
    ),
)
def run(
    context: typer.Context,
    assembly: Annotated[
        Path,
        typer.Option(
            "--assembly",
            "-i",
            help="Assembled nucleotide FASTA to process.",
            exists=True,
            dir_okay=False,
            readable=True,
            rich_help_panel="Inputs",
        ),
    ],
    database: Annotated[
        Path,
        typer.Option(
            "--database",
            "-d",
            help="Database directory produced by 'petra build-db'.",
            exists=True,
            file_okay=False,
            readable=True,
            callback=_validate_database_directory,
            rich_help_panel="Inputs",
        ),
    ],
    long_reads: Annotated[
        list[Path],
        typer.Option(
            "--long-reads",
            "-r",
            help="Long-read sequence file to map; repeat for multiple files.",
            exists=True,
            dir_okay=False,
            readable=True,
            rich_help_panel="Inputs",
        ),
    ],
    technology: Annotated[
        SequencingTechnology,
        typer.Option(
            "--technology",
            help="Long-read sequencing technology used to select the Minimap2 preset.",
            rich_help_panel="Long-read mapping",
            show_default=True,
        ),
    ] = SequencingTechnology.PACBIO,
    no_bam_cleaning: Annotated[
        bool,
        typer.Option(
            "--no-bam-cleaning",
            help="Keep the sorted BAM after median coverage has been calculated.",
            rich_help_panel="Long-read mapping",
        ),
    ] = False,
    max_low_complexity: Annotated[
        int,
        typer.Option(
            "--max-low-complexity",
            min=0,
            help=(
                "Maximum number of low-complexity amino-acid residues allowed "
                "per predicted protein."
            ),
            rich_help_panel="Gene prediction",
            show_default=True,
        ),
    ] = 100,
    max_clusters_per_sequence: Annotated[
        int,
        typer.Option(
            "--max-clusters",
            min=1,
            help=(
                "Maximum number of protein clusters retained per input "
                "assembly sequence."
            ),
            rich_help_panel="Gene prediction",
            show_default=True,
        ),
    ] = 1000,
    max_target_occurrences: Annotated[
        int,
        typer.Option(
            "--max-target-occurrences",
            min=0,
            help=(
                "Maximum number of alignments allowed per target sequence. "
                "Targets above this count are removed."
            ),
            rich_help_panel="Protein alignment",
            show_default=True,
        ),
    ] = 50,
    alpha: Annotated[
        float,
        typer.Option(
            "--alpha",
            min=0.0,
            max=1.0,
            callback=_validate_alpha,
            help=(
                "Minimum share of a sequence's aligned clusters that must support "
                "a taxon for it to be assigned. The share must hold at the lower "
                "bound of its 95% confidence interval, so sequences with few "
                "aligned clusters stay unassigned."
            ),
            rich_help_panel="Taxonomic assignment",
            show_default=True,
        ),
    ] = 0.5,
    rrna: Annotated[
        bool,
        typer.Option(
            "--rrna",
            help=(
                "Detect rRNA genes with cmscan and classify SSU and LSU genes "
                "with blastn."
            ),
            rich_help_panel="rRNA detection",
        ),
    ] = False,
    rrna_min_alignment_length: Annotated[
        int,
        typer.Option(
            "--rrna-min-alignment-length",
            min=1,
            help="Minimum BLAST alignment length for a hit to be reported.",
            rich_help_panel="rRNA detection",
            show_default=True,
        ),
    ] = 200,
    host_taxon: Annotated[
        str | None,
        typer.Option(
            "--host-taxon",
            help=(
                "Expected host taxon, as an NCBI scientific name or TaxID. rRNA "
                "hits whose lineage does not contain it are flagged non_host."
            ),
            rich_help_panel="rRNA detection",
        ),
    ] = None,
    trna: Annotated[
        bool,
        typer.Option(
            "--trna",
            help="Annotate tRNA genes with tRNAscan-SE.",
            rich_help_panel="tRNA detection",
        ),
    ] = False,
    summary_font_size: Annotated[
        float,
        typer.Option(
            "--summary-font-size",
            min=1.0,
            help="Base font size of the summary graph in points.",
            rich_help_panel="Summary graph",
            show_default=True,
        ),
    ] = 10.0,
    summary_colors: Annotated[
        Path | None,
        typer.Option(
            "--summary-colors",
            exists=True,
            dir_okay=False,
            readable=True,
            help=(
                "Two-column TSV of taxonomic rank and color. Sequences are colored "
                "with the deepest assigned rank found in the file instead of their "
                "superkingdom; assigned sequences without a match are colored as "
                "'other'."
            ),
            rich_help_panel="Summary graph",
        ),
    ] = None,
    summary_disable_panel: Annotated[
        list[SummaryPanel] | None,
        typer.Option(
            "--summary-disable-panel",
            help="Summary graph panel to leave out; repeat to drop several.",
            rich_help_panel="Summary graph",
        ),
    ] = None,
    output: Annotated[
        Path,
        typer.Option(
            "--output",
            "-o",
            help="Required output directory and Snakemake working directory.",
            file_okay=False,
            rich_help_panel="Execution",
        ),
    ] = Path("petra_output"),
    profile: Annotated[
        str,
        typer.Option(
            "--profile",
            help="Packaged profile name or external workflow-profile path.",
            rich_help_panel="Execution",
            show_default=True,
        ),
    ] = "default",
    executor: Annotated[
        str | None,
        typer.Option(
            "--executor",
            help="Override the selected profile's executor.",
            rich_help_panel="Execution",
        ),
    ] = "slurm",
) -> None:
    """Run the complete workflow unless Snakemake targets are passed after '--'."""
    if host_taxon is not None and not rrna:
        raise typer.BadParameter(
            "--host-taxon requires --rrna.",
            param_hint="--host-taxon",
        )

    if rrna:
        rna_database_filenames = (
            COVARIANCE_MODEL_FILENAME,
            CLAN_INFO_FILENAME,
            SSU_DATABASE_FILENAME,
            LSU_DATABASE_FILENAME,
        )
        missing = [
            filename
            for filename in rna_database_filenames
            if not (database / filename).is_file()
        ]
        if missing:
            raise typer.BadParameter(
                f"RNA database files not found in {database}: {', '.join(missing)}. "
                "Build them with "
                f"'petra build-db --output {database} --build-rna-dbs "
                "--skip-diamond-db'.",
                param_hint="--database",
            )

    try:
        # taxonkit resolves the host taxon, so PATH is checked first to report a
        # missing tool as such instead of as a resolution failure.
        _require_executables(rrna=rrna, trna=trna)
        canonical_host_taxon = (
            resolve_host_taxon(host_taxon) if host_taxon is not None else None
        )
        workflow_config = {
            "assembly": assembly,
            "diamond_db": database / DIAMOND_DATABASE_FILENAME,
            "long_reads": long_reads,
            "sequencing_technology": technology,
            "no_bam_cleaning": no_bam_cleaning,
            "max_low_complexity": max_low_complexity,
            "max_clusters_per_sequence": max_clusters_per_sequence,
            "max_target_occurrences": max_target_occurrences,
            "alpha": alpha,
            "summary_font_size": summary_font_size,
            "summary_colors": summary_colors,
            "summary_disabled_panels": [
                panel.value for panel in summary_disable_panel or ()
            ],
            "rrna": rrna,
            "trna": trna,
            "rrna_chunks": _RRNA_CHUNKS,
            "rrna_min_alignment_length": rrna_min_alignment_length,
            "host_taxon": canonical_host_taxon,
            "covariance_model": database / COVARIANCE_MODEL_FILENAME,
            "clan_info": database / CLAN_INFO_FILENAME,
            "ssu_database": database / SSU_DATABASE_FILENAME,
            "lsu_database": database / LSU_DATABASE_FILENAME,
        }
        exit_code = launch_workflow(
            distribution_name="petra",
            output=output,
            workflow_config=workflow_config,
            profile=profile,
            executor=executor,
            passthrough=context.args,
            console=console,
        )
    except (HostTaxonError, LauncherError) as error:
        typer.echo(f"Error: {error}", err=True)
        raise typer.Exit(code=1) from error

    raise typer.Exit(code=exit_code)


@app.command(
    "build-db",
    context_settings=_HELP_CONTEXT,
    epilog=(
        "The DIAMOND database usually takes around 11 hours and 20GB of RAM to "
        "complete. The rRNA databases take about an hour."
    ),
)
def build_db(
    output: Annotated[
        Path,
        typer.Option(
            "--output",
            "-o",
            help="Directory in which to create the databases.",
            file_okay=False,
            rich_help_panel="Output",
        ),
    ],
    fasta: Annotated[
        Path | None,
        typer.Option(
            "--fasta",
            help="Protein FASTA input. Downloads NCBI NR when omitted.",
            dir_okay=False,
            rich_help_panel="Inputs",
        ),
    ] = None,
    taxon_map: Annotated[
        Path | None,
        typer.Option(
            "--taxon-map",
            help=(
                "NCBI protein accession-to-taxon map. Downloads the full map "
                "when omitted."
            ),
            dir_okay=False,
            rich_help_panel="Inputs",
        ),
    ] = None,
    build_rna_dbs: Annotated[
        bool,
        typer.Option(
            "--build-rna-dbs",
            help=(
                "Also build the rRNA covariance model and the SSU and LSU BLAST "
                "databases next to petra.dmnd. Requires makeblastdb, cmpress and "
                "taxonkit on PATH."
            ),
            rich_help_panel="Databases",
        ),
    ] = False,
    skip_diamond_db: Annotated[
        bool,
        typer.Option(
            "--skip-diamond-db",
            help=(
                "Leave petra.dmnd alone, to build the rRNA databases in a second "
                "step without repeating the DIAMOND build."
            ),
            rich_help_panel="Databases",
        ),
    ] = False,
) -> None:
    """Build the taxonomy-aware DIAMOND database and, optionally, the rRNA ones."""
    if skip_diamond_db and not build_rna_dbs:
        raise typer.BadParameter(
            "--skip-diamond-db leaves nothing to build. Add --build-rna-dbs.",
            param_hint="--skip-diamond-db",
        )
    if skip_diamond_db and (fasta is not None or taxon_map is not None):
        raise typer.BadParameter(
            "--fasta and --taxon-map only feed the DIAMOND database, which "
            "--skip-diamond-db disables.",
            param_hint="--skip-diamond-db",
        )

    try:
        if not skip_diamond_db:
            build_diamond_database(
                fasta=fasta,
                taxon_map=taxon_map,
                output=output,
            )
        if build_rna_dbs:
            build_rna_databases(output=output)
    except DatabaseBuildError as error:
        typer.echo(f"Error: {error}", err=True)
        raise typer.Exit(code=1) from error


@app.command("profiles", context_settings=_HELP_CONTEXT)
def profiles() -> None:
    """List packaged workflow profiles and common Slurm defaults."""
    try:
        summaries = profile_summaries()
    except LauncherError as error:
        typer.echo(f"Error: {error}", err=True)
        raise typer.Exit(code=1) from error

    if not summaries:
        typer.echo("No packaged workflow profiles found.")
        return

    table = Table(show_header=True, header_style="bold")
    table.add_column("PROFILE")
    table.add_column("JOBS")
    table.add_column("CORES")
    table.add_column("PARTITION")
    table.add_column("QOS")
    table.add_column("RUNTIME")
    table.add_column("MEM/CPU")
    table.add_column("PATH")
    for summary in summaries:
        table.add_row(
            summary.name,
            summary.jobs,
            summary.cores,
            summary.partition,
            summary.qos,
            summary.runtime,
            summary.memory_per_cpu,
            summary.path,
        )
    console.print(table)
