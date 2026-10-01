from __future__ import annotations

import logging
import shutil
from enum import Enum
from importlib.metadata import PackageNotFoundError, version
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
    TAXONOMY_DIRECTORY_NAME,
    DatabaseBuildError,
    build_diamond_database,
    build_rna_databases,
    install_taxonomy_database,
    validate_taxonomy_directory,
)
from .launcher import LauncherError, launch_workflow, profile_summaries
from .taxonomy import HostTaxonError, resolve_host_phylum

_HELP_CONTEXT = {"help_option_names": ["-h", "--help"]}
_RUN_CONTEXT = {**_HELP_CONTEXT, "allow_extra_args": True}
_RRNA_CHUNKS = 4
_THEMES_DIRECTORY = Path(__file__).resolve().parent / "themes"


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


SummaryTheme = Enum(
    "SummaryTheme",
    {path.stem: path.stem for path in sorted(_THEMES_DIRECTORY.glob("*.tsv"))},
    type=str,
)


def _require_executables(
    *, rrna: bool, trna: bool, skip_taxonomic_assignment: bool
) -> None:
    required = ["pyrodigal", "whokaryote.py"]
    if not skip_taxonomic_assignment:
        required.extend(
            (
                "diamond",
                "taxonkit",
                "cd-hit",
                "minimap2",
                "samtools",
                "segmasker",
                "bedtools",
            )
        )
    if rrna:
        required.extend(("cmscan", "blastn"))
        if skip_taxonomic_assignment:
            required.append("taxonkit")
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
    if not 0 < value < 1:
        raise typer.BadParameter(
            "The alpha-majority threshold must be greater than 0 and less than 1."
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
        Path | None,
        typer.Option(
            "--database",
            "-d",
            help=(
                "Database directory produced by 'petra build-db'. Required for "
                "taxonomic assignment and rRNA analysis. Can also be set with "
                "PETRA_DB_DIR; --database takes precedence."
            ),
            envvar="PETRA_DB_DIR",
            show_envvar=True,
            exists=True,
            file_okay=False,
            readable=True,
            rich_help_panel="Inputs",
        ),
    ] = None,
    long_reads: Annotated[
        list[Path] | None,
        typer.Option(
            "--long-reads",
            "-r",
            help=(
                "Long-read sequence file to map; repeat for multiple files. "
                "Required unless taxonomic assignment is skipped."
            ),
            exists=True,
            dir_okay=False,
            readable=True,
            rich_help_panel="Inputs",
        ),
    ] = None,
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
    max_unigenes_per_sequence: Annotated[
        int,
        typer.Option(
            "--max-unigenes",
            min=1,
            help=(
                "Maximum number of unigenes retained per input "
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
            callback=_validate_alpha,
            help=(
                "Share of a sequence's aligned unigenes that a taxon must "
                "exceed, strictly between 0 and 1; the default requires a "
                "strict majority."
            ),
            rich_help_panel="Taxonomic assignment",
            show_default=True,
        ),
    ] = 0.5,
    min_aligned_unigenes: Annotated[
        int,
        typer.Option(
            "--min-aligned-unigenes",
            min=3,
            help=(
                "Minimum aligned unigenes required per sequence at or below "
                "--unigene-support-length. Longer sequences require more support."
            ),
            rich_help_panel="Taxonomic assignment",
            show_default=True,
        ),
    ] = 3,
    unigene_support_length: Annotated[
        int,
        typer.Option(
            "--unigene-support-length",
            min=1,
            help="Sequence length in bp above which required unigene support grows.",
            rich_help_panel="Taxonomic assignment",
            show_default=True,
        ),
    ] = 50_000,
    unigenes_per_decade: Annotated[
        int,
        typer.Option(
            "--unigenes-per-decade",
            min=1,
            help=(
                "Additional aligned unigenes required per tenfold length increase "
                "above --unigene-support-length."
            ),
            rich_help_panel="Taxonomic assignment",
            show_default=True,
        ),
    ] = 3,
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
                "Expected host taxon, as an NCBI scientific name or TaxID. PETRA "
                "resolves its phylum and flags rRNA hits outside that phylum. "
                "If no phylum is available, reports non_host=- and writes "
                "rna/predictions.rrna.warnings.txt."
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
    skip_taxonomic_assignment: Annotated[
        bool,
        typer.Option(
            "--skip-taxonomic-assignment",
            help=(
                "Skip the protein-based taxonomic assignment and run only "
                "Whokaryote plus the enabled rRNA and tRNA analyses."
            ),
            rich_help_panel="Execution",
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
                "with the deepest assigned rank found in the file; assigned "
                "sequences without a match are colored as 'other'. Cannot be "
                "combined with --summary-theme."
            ),
            rich_help_panel="Summary graph",
        ),
    ] = None,
    summary_theme: Annotated[
        SummaryTheme | None,
        typer.Option(
            "--summary-theme",
            help=(
                "Packaged color file of the summary graph. 'superkingdom' is used "
                "unless this or --summary-colors is given."
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
    """Run the workflow unless Snakemake targets are passed after '--'."""
    if host_taxon is not None and not rrna:
        raise typer.BadParameter(
            "--host-taxon requires --rrna.",
            param_hint="--host-taxon",
        )
    if summary_colors is not None and summary_theme is not None:
        raise typer.BadParameter(
            "--summary-colors cannot be combined with --summary-theme.",
            param_hint="--summary-theme",
        )
    if database is None and (not skip_taxonomic_assignment or rrna):
        raise typer.BadParameter(
            "--database is required for taxonomic assignment and rRNA analysis.",
            param_hint="--database",
        )
    if not skip_taxonomic_assignment and not long_reads:
        raise typer.BadParameter(
            "--long-reads is required unless --skip-taxonomic-assignment is used.",
            param_hint="--long-reads",
        )
    if (
        database is not None
        and not skip_taxonomic_assignment
        and not (database / DIAMOND_DATABASE_FILENAME).is_file()
    ):
        raise typer.BadParameter(
            f"{DIAMOND_DATABASE_FILENAME} not found in {database}. Build it with "
            f"'petra build-db --output {database}'.",
            param_hint="--database",
        )

    if rrna:
        assert database is not None
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
    taxonomy = Path()
    if not skip_taxonomic_assignment or rrna:
        assert database is not None
        try:
            taxonomy = validate_taxonomy_directory(database / TAXONOMY_DIRECTORY_NAME)
        except DatabaseBuildError as error:
            raise typer.BadParameter(
                f"{error} Build it with 'petra build-db --output {database}'.",
                param_hint="--database",
            ) from error
    try:
        # taxonkit resolves the host phylum, so PATH is checked first to report a
        # missing tool as such instead of as a resolution failure.
        _require_executables(
            rrna=rrna,
            trna=trna,
            skip_taxonomic_assignment=skip_taxonomic_assignment,
        )
        host_phylum = (
            resolve_host_phylum(host_taxon, taxonomy=taxonomy)
            if host_taxon is not None
            else None
        )
        database_path = database if database is not None else Path()
        workflow_config = {
            "assembly": assembly,
            "diamond_db": database_path / DIAMOND_DATABASE_FILENAME,
            "long_reads": long_reads or [],
            "sequencing_technology": technology,
            "no_bam_cleaning": no_bam_cleaning,
            "max_low_complexity": max_low_complexity,
            "max_unigenes_per_sequence": max_unigenes_per_sequence,
            "max_target_occurrences": max_target_occurrences,
            "alpha": alpha,
            "min_aligned_unigenes": min_aligned_unigenes,
            "unigene_support_length": unigene_support_length,
            "unigenes_per_decade": unigenes_per_decade,
            "summary_font_size": summary_font_size,
            "summary_colors": (
                _THEMES_DIRECTORY / f"{summary_theme.value}.tsv"
                if summary_theme is not None
                else summary_colors
            ),
            "summary_disabled_panels": [
                panel.value for panel in summary_disable_panel or ()
            ],
            "rrna": rrna,
            "trna": trna,
            "skip_taxonomic_assignment": skip_taxonomic_assignment,
            "rrna_chunks": _RRNA_CHUNKS,
            "rrna_min_alignment_length": rrna_min_alignment_length,
            "host_taxon": host_taxon,
            "host_phylum": host_phylum,
            "covariance_model": database_path / COVARIANCE_MODEL_FILENAME,
            "clan_info": database_path / CLAN_INFO_FILENAME,
            "ssu_database": database_path / SSU_DATABASE_FILENAME,
            "lsu_database": database_path / LSU_DATABASE_FILENAME,
            "taxonomy": taxonomy,
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
    taxdump: Annotated[
        Path | None,
        typer.Option(
            "--taxdump",
            help=(
                "Local NCBI taxdump archive. Downloads the current archive when "
                "omitted."
            ),
            exists=True,
            dir_okay=False,
            readable=True,
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
    """Build PETRA databases with their matching NCBI taxonomy data."""
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
        taxonomy = install_taxonomy_database(output=output, taxdump=taxdump)
        if not skip_diamond_db:
            build_diamond_database(
                fasta=fasta,
                taxon_map=taxon_map,
                output=output,
            )
        if build_rna_dbs:
            build_rna_databases(output=output, taxonomy=taxonomy)
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
