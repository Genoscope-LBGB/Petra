from __future__ import annotations

import logging
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.logging import RichHandler
from rich.table import Table

from .database import DatabaseBuildError, build_diamond_database
from .launcher import LauncherError, launch_workflow, profile_summaries

_HELP_CONTEXT = {"help_option_names": ["-h", "--help"]}
_RUN_CONTEXT = {**_HELP_CONTEXT, "allow_extra_args": True}


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
        "petra run --output run-1 -- --targets hello"
    ),
)
def run(
    context: typer.Context,
    output: Annotated[
        Path,
        typer.Option(
            "--output",
            "-o",
            help="Required output directory and Snakemake working directory.",
            file_okay=False,
            rich_help_panel="Execution",
        ),
    ],
    message: Annotated[
        str,
        typer.Option(
            "--message",
            help="Text written by the starter workflow.",
            rich_help_panel="Workflow options",
            show_default=True,
        ),
    ] = "Hello from petra",
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
    profile_source = context.get_parameter_source("profile")

    workflow_config = {
        "message": message,
    }
    try:
        exit_code = launch_workflow(
            distribution_name="petra",
            output=output,
            workflow_config=workflow_config,
            profile=profile,
            executor=executor,
            passthrough=context.args,
            console=console,
        )
    except LauncherError as error:
        typer.echo(f"Error: {error}", err=True)
        raise typer.Exit(code=1) from error

    raise typer.Exit(code=exit_code)


@app.command("build-db", context_settings=_HELP_CONTEXT)
def build_db(
    output: Annotated[
        Path,
        typer.Option(
            "--output",
            "-o",
            help="Directory in which to create petra.dmnd.",
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
) -> None:
    """Build a taxonomy-aware DIAMOND database."""
    try:
        build_diamond_database(
            fasta=fasta,
            taxon_map=taxon_map,
            output=output,
        )
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
