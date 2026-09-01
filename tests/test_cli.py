from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from typer.testing import CliRunner

from petra.cli import SequencingTechnology, app
from petra.database import DatabaseBuildError


class RunCliTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.database = self.directory / "database"
        self.database.mkdir()
        (self.database / "petra.dmnd").touch()
        self.assembly = self.directory / "assembly.fasta"
        self.reads = self.directory / "reads.fastq"
        for path in (self.assembly, self.reads):
            path.touch()

    def test_accepts_repeated_long_reads_and_ont_technology(self) -> None:
        second_reads = self.directory / "second.fastq"
        second_reads.touch()

        with (
            patch("petra.cli._require_executables"),
            patch("petra.cli.launch_workflow", return_value=0) as launch,
        ):
            result = CliRunner().invoke(
                app,
                [
                    "run",
                    "--assembly",
                    str(self.assembly),
                    "--database",
                    str(self.database),
                    "--long-reads",
                    str(self.reads),
                    "--long-reads",
                    str(second_reads),
                    "--technology",
                    "ont",
                ],
            )

        self.assertEqual(result.exit_code, 0, result.output)
        workflow_config = launch.call_args.kwargs["workflow_config"]
        self.assertEqual(workflow_config["long_reads"], [self.reads, second_reads])
        self.assertEqual(
            workflow_config["sequencing_technology"], SequencingTechnology.ONT
        )

    def test_resolves_the_diamond_database_inside_the_directory(self) -> None:
        with (
            patch("petra.cli._require_executables"),
            patch("petra.cli.launch_workflow", return_value=0) as launch,
        ):
            result = CliRunner().invoke(
                app,
                [
                    "run",
                    "-i",
                    str(self.assembly),
                    "-d",
                    str(self.database),
                    "-r",
                    str(self.reads),
                ],
            )

        self.assertEqual(result.exit_code, 0, result.output)
        workflow_config = launch.call_args.kwargs["workflow_config"]
        self.assertEqual(
            workflow_config["diamond_db"], self.database / "petra.dmnd"
        )
        self.assertEqual(
            workflow_config["sequencing_technology"], SequencingTechnology.PACBIO
        )

    def test_rejects_a_directory_without_the_diamond_database(self) -> None:
        empty = self.directory / "empty"
        empty.mkdir()

        with (
            patch("petra.cli._require_executables"),
            patch("petra.cli.launch_workflow") as launch,
        ):
            result = CliRunner().invoke(
                app,
                [
                    "run",
                    "--assembly",
                    str(self.assembly),
                    "--database",
                    str(empty),
                    "--long-reads",
                    str(self.reads),
                ],
            )

        self.assertEqual(result.exit_code, 2)
        self.assertIn(
            f"petra.dmnd not found in {empty}. Build it with 'petra build-db "
            f"--output {empty}'.",
            _flatten(result.output),
        )
        launch.assert_not_called()

    def test_rejects_a_diamond_database_file_as_the_database(self) -> None:
        with (
            patch("petra.cli._require_executables"),
            patch("petra.cli.launch_workflow") as launch,
        ):
            result = CliRunner().invoke(
                app,
                [
                    "run",
                    "--assembly",
                    str(self.assembly),
                    "--database",
                    str(self.database / "petra.dmnd"),
                    "--long-reads",
                    str(self.reads),
                ],
            )

        self.assertEqual(result.exit_code, 2)
        launch.assert_not_called()

    def test_rejects_missing_executables_before_launch(self) -> None:
        def which_without_bedtools(executable: str) -> str | None:
            return None if executable == "bedtools" else f"/usr/bin/{executable}"

        with (
            patch("petra.cli.shutil.which", side_effect=which_without_bedtools),
            patch("petra.cli.launch_workflow") as launch,
        ):
            result = CliRunner().invoke(
                app,
                [
                    "run",
                    "--assembly",
                    str(self.assembly),
                    "--database",
                    str(self.database),
                    "--long-reads",
                    str(self.reads),
                ],
            )

        self.assertEqual(result.exit_code, 1)
        self.assertIn(
            "Required executables not found on PATH: bedtools.",
            result.output,
        )
        launch.assert_not_called()


def _flatten(output: str) -> str:
    """Collapse the box drawing and line wrapping Rich applies to usage errors."""
    return " ".join(output.replace("│", " ").split())


class BuildDbCliTest(unittest.TestCase):
    def test_builds_only_the_diamond_database_by_default(self) -> None:
        with (
            patch("petra.cli.build_diamond_database") as diamond,
            patch("petra.cli.build_rna_databases") as rna,
        ):
            result = CliRunner().invoke(app, ["build-db", "--output", "database"])

        self.assertEqual(result.exit_code, 0, result.output)
        diamond.assert_called_once()
        rna.assert_not_called()

    def test_builds_both_databases_in_the_same_directory(self) -> None:
        with (
            patch("petra.cli.build_diamond_database") as diamond,
            patch("petra.cli.build_rna_databases") as rna,
        ):
            result = CliRunner().invoke(
                app, ["build-db", "--output", "database", "--build-rna-dbs"]
            )

        self.assertEqual(result.exit_code, 0, result.output)
        self.assertEqual(
            diamond.call_args.kwargs["output"], rna.call_args.kwargs["output"]
        )

    def test_skips_the_diamond_database_when_asked(self) -> None:
        with (
            patch("petra.cli.build_diamond_database") as diamond,
            patch("petra.cli.build_rna_databases") as rna,
        ):
            result = CliRunner().invoke(
                app,
                [
                    "build-db",
                    "--output",
                    "database",
                    "--build-rna-dbs",
                    "--skip-diamond-db",
                ],
            )

        self.assertEqual(result.exit_code, 0, result.output)
        diamond.assert_not_called()
        rna.assert_called_once_with(output=Path("database"))

    def test_rejects_skipping_every_database(self) -> None:
        with (
            patch("petra.cli.build_diamond_database") as diamond,
            patch("petra.cli.build_rna_databases") as rna,
        ):
            result = CliRunner().invoke(
                app, ["build-db", "--output", "database", "--skip-diamond-db"]
            )

        self.assertEqual(result.exit_code, 2)
        self.assertIn(
            "--skip-diamond-db leaves nothing to build. Add --build-rna-dbs.",
            _flatten(result.output),
        )
        diamond.assert_not_called()
        rna.assert_not_called()

    def test_rejects_diamond_inputs_when_the_diamond_build_is_skipped(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            fasta = Path(temporary_directory) / "proteins.faa"
            fasta.touch()

            with (
                patch("petra.cli.build_diamond_database") as diamond,
                patch("petra.cli.build_rna_databases") as rna,
            ):
                result = CliRunner().invoke(
                    app,
                    [
                        "build-db",
                        "--output",
                        "database",
                        "--build-rna-dbs",
                        "--skip-diamond-db",
                        "--fasta",
                        str(fasta),
                    ],
                )

        self.assertEqual(result.exit_code, 2)
        self.assertIn(
            "--fasta and --taxon-map only feed the DIAMOND database",
            _flatten(result.output),
        )
        diamond.assert_not_called()
        rna.assert_not_called()

    def test_reports_a_build_failure(self) -> None:
        with (
            patch("petra.cli.build_diamond_database"),
            patch(
                "petra.cli.build_rna_databases",
                side_effect=DatabaseBuildError("cmpress failed with exit code 1."),
            ),
        ):
            result = CliRunner().invoke(
                app, ["build-db", "--output", "database", "--build-rna-dbs"]
            )

        self.assertEqual(result.exit_code, 1)
        self.assertIn("cmpress failed with exit code 1.", result.output)


if __name__ == "__main__":
    unittest.main()
