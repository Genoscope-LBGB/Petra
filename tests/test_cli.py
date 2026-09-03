from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from typer.testing import CliRunner

from petra.cli import SequencingTechnology, app
from petra.database import DatabaseBuildError
from petra.taxonomy import HostTaxonError


class RunCliTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.environment = patch.dict("os.environ", {"PETRA_DB_DIR": ""})
        self.environment.start()
        self.addCleanup(self.environment.stop)
        self.directory = Path(self.temporary.name)
        self.database = self.directory / "database"
        self.database.mkdir()
        (self.database / "petra.dmnd").touch()
        for filename in (
            "petra.rrna.cm",
            "petra.rrna.clanin",
            "petra.ssu.fasta",
            "petra.lsu.fasta",
        ):
            (self.database / filename).touch()
        self.assembly = self.directory / "assembly.fasta"
        self.reads = self.directory / "reads.fastq"
        for path in (self.assembly, self.reads):
            path.touch()

    def run_arguments(self, *additional: str) -> list[str]:
        return [
            "run",
            "--assembly",
            str(self.assembly),
            "--database",
            str(self.database),
            "--long-reads",
            str(self.reads),
            *additional,
        ]

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
        self.assertEqual(workflow_config["diamond_db"], self.database / "petra.dmnd")
        self.assertEqual(
            workflow_config["sequencing_technology"], SequencingTechnology.PACBIO
        )

    def test_uses_database_directory_from_environment(self) -> None:
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
                    "--long-reads",
                    str(self.reads),
                ],
                env={"PETRA_DB_DIR": str(self.database)},
            )

        self.assertEqual(result.exit_code, 0, result.output)
        workflow_config = launch.call_args.kwargs["workflow_config"]
        self.assertEqual(workflow_config["diamond_db"], self.database / "petra.dmnd")

    def test_database_option_overrides_environment(self) -> None:
        environment_database = self.directory / "environment-database"
        environment_database.mkdir()
        (environment_database / "petra.dmnd").touch()

        with (
            patch("petra.cli._require_executables"),
            patch("petra.cli.launch_workflow", return_value=0) as launch,
        ):
            result = CliRunner().invoke(
                app,
                self.run_arguments(),
                env={"PETRA_DB_DIR": str(environment_database)},
            )

        self.assertEqual(result.exit_code, 0, result.output)
        workflow_config = launch.call_args.kwargs["workflow_config"]
        self.assertEqual(workflow_config["diamond_db"], self.database / "petra.dmnd")

    def test_checks_database_files_for_environment_directory(self) -> None:
        empty = self.directory / "empty-environment-database"
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
                    "--long-reads",
                    str(self.reads),
                ],
                env={"PETRA_DB_DIR": str(empty)},
            )

        self.assertEqual(result.exit_code, 2)
        self.assertIn(
            f"petra.dmnd not found in {empty}. Build it with 'petra build-db "
            f"--output {empty}'.",
            _flatten(result.output),
        )
        launch.assert_not_called()

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

    def test_rrna_populates_workflow_config(self) -> None:
        with (
            patch("petra.cli._require_executables"),
            patch("petra.cli.launch_workflow", return_value=0) as launch,
        ):
            result = CliRunner().invoke(app, self.run_arguments("--rrna"))

        self.assertEqual(result.exit_code, 0, result.output)
        workflow_config = launch.call_args.kwargs["workflow_config"]
        self.assertIs(workflow_config["rrna"], True)
        self.assertIs(workflow_config["trna"], False)
        self.assertEqual(workflow_config["rrna_chunks"], 4)
        self.assertEqual(workflow_config["rrna_min_alignment_length"], 200)
        self.assertIsNone(workflow_config["host_phylum"])
        self.assertIs(workflow_config["skip_taxonomic_assignment"], False)
        self.assertEqual(
            workflow_config["covariance_model"],
            self.database / "petra.rrna.cm",
        )
        self.assertEqual(
            workflow_config["clan_info"],
            self.database / "petra.rrna.clanin",
        )
        self.assertEqual(
            workflow_config["ssu_database"],
            self.database / "petra.ssu.fasta",
        )
        self.assertEqual(
            workflow_config["lsu_database"],
            self.database / "petra.lsu.fasta",
        )

    def test_skips_taxonomic_assignment_for_rrna_without_reads_or_diamond(self) -> None:
        (self.database / "petra.dmnd").unlink()

        def which_with_only_rrna_tools(executable: str) -> str | None:
            return (
                f"/usr/bin/{executable}"
                if executable in {"cmscan", "blastn", "taxonkit"}
                else None
            )

        with (
            patch("petra.cli.shutil.which", side_effect=which_with_only_rrna_tools),
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
                    "--rrna",
                    "--skip-taxonomic-assignment",
                ],
            )

        self.assertEqual(result.exit_code, 0, result.output)
        workflow_config = launch.call_args.kwargs["workflow_config"]
        self.assertIs(workflow_config["skip_taxonomic_assignment"], True)
        self.assertIs(workflow_config["rrna"], True)
        self.assertIs(workflow_config["trna"], False)
        self.assertEqual(workflow_config["long_reads"], [])

    def test_skips_taxonomic_assignment_for_trna_without_reads_or_database(
        self,
    ) -> None:
        def which_with_only_trna_tools(executable: str) -> str | None:
            return f"/usr/bin/{executable}" if executable == "tRNAscan-SE" else None

        with (
            patch("petra.cli.shutil.which", side_effect=which_with_only_trna_tools),
            patch("petra.cli.launch_workflow", return_value=0) as launch,
        ):
            result = CliRunner().invoke(
                app,
                [
                    "run",
                    "--assembly",
                    str(self.assembly),
                    "--trna",
                    "--skip-taxonomic-assignment",
                ],
            )

        self.assertEqual(result.exit_code, 0, result.output)
        workflow_config = launch.call_args.kwargs["workflow_config"]
        self.assertIs(workflow_config["skip_taxonomic_assignment"], True)
        self.assertIs(workflow_config["rrna"], False)
        self.assertIs(workflow_config["trna"], True)
        self.assertEqual(workflow_config["long_reads"], [])

    def test_rejects_skip_taxonomic_assignment_without_rna_analysis(self) -> None:
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
                    "--skip-taxonomic-assignment",
                ],
            )

        self.assertEqual(result.exit_code, 2)
        self.assertIn(
            "--skip-taxonomic-assignment requires --rrna, --trna, or both.",
            _flatten(result.output),
        )
        launch.assert_not_called()

    def test_rejects_rrna_when_an_rna_database_file_is_missing(self) -> None:
        (self.database / "petra.rrna.clanin").unlink()

        with (
            patch("petra.cli._require_executables"),
            patch("petra.cli.launch_workflow") as launch,
        ):
            result = CliRunner().invoke(app, self.run_arguments("--rrna"))

        self.assertEqual(result.exit_code, 2)
        output = _flatten(result.output)
        self.assertIn("petra.rrna.clanin", output)
        self.assertIn(
            f"petra build-db --output {self.database} --build-rna-dbs "
            "--skip-diamond-db",
            output,
        )
        launch.assert_not_called()

    def test_rejects_host_taxon_without_rrna(self) -> None:
        with (
            patch("petra.cli._require_executables"),
            patch("petra.cli.launch_workflow") as launch,
        ):
            result = CliRunner().invoke(
                app,
                self.run_arguments("--host-taxon", "9606"),
            )

        self.assertEqual(result.exit_code, 2)
        self.assertIn("--host-taxon requires --rrna.", _flatten(result.output))
        launch.assert_not_called()

    def test_resolves_host_taxon_phylum_into_workflow_config(self) -> None:
        with (
            patch(
                "petra.cli.resolve_host_phylum",
                return_value="Chordata",
            ) as resolve,
            patch("petra.cli._require_executables"),
            patch("petra.cli.launch_workflow", return_value=0) as launch,
        ):
            result = CliRunner().invoke(
                app,
                self.run_arguments("--rrna", "--host-taxon", "9606"),
            )

        self.assertEqual(result.exit_code, 0, result.output)
        resolve.assert_called_once_with("9606")
        self.assertEqual(
            launch.call_args.kwargs["workflow_config"]["host_phylum"],
            "Chordata",
        )

    def test_reports_host_taxon_resolution_failure(self) -> None:
        with (
            patch(
                "petra.cli.resolve_host_phylum",
                side_effect=HostTaxonError("TaxID 123 is not usable."),
            ),
            patch("petra.cli._require_executables"),
            patch("petra.cli.launch_workflow") as launch,
        ):
            result = CliRunner().invoke(
                app,
                self.run_arguments("--rrna", "--host-taxon", "123"),
            )

        self.assertEqual(result.exit_code, 1)
        self.assertIn("Error: TaxID 123 is not usable.", result.output)
        launch.assert_not_called()

    def test_trna_requires_trnascan_se(self) -> None:
        def which_without_trnascan(executable: str) -> str | None:
            return None if executable == "tRNAscan-SE" else f"/usr/bin/{executable}"

        with (
            patch("petra.cli.shutil.which", side_effect=which_without_trnascan),
            patch("petra.cli.launch_workflow") as launch,
        ):
            result = CliRunner().invoke(app, self.run_arguments("--trna"))

        self.assertEqual(result.exit_code, 1)
        self.assertIn(
            "Required executables not found on PATH: tRNAscan-SE.",
            result.output,
        )
        launch.assert_not_called()

    def test_rrna_requires_cmscan_and_blastn(self) -> None:
        def which_without_rrna_tools(executable: str) -> str | None:
            return (
                None if executable in {"cmscan", "blastn"} else f"/usr/bin/{executable}"
            )

        with (
            patch("petra.cli.shutil.which", side_effect=which_without_rrna_tools),
            patch("petra.cli.launch_workflow") as launch,
        ):
            result = CliRunner().invoke(app, self.run_arguments("--rrna"))

        self.assertEqual(result.exit_code, 1)
        self.assertIn(
            "Required executables not found on PATH: cmscan, blastn.",
            result.output,
        )
        launch.assert_not_called()

    def test_plain_run_does_not_require_rna_tools(self) -> None:
        def which_without_rna_tools(executable: str) -> str | None:
            return (
                None
                if executable in {"cmscan", "blastn", "tRNAscan-SE"}
                else f"/usr/bin/{executable}"
            )

        with (
            patch("petra.cli.shutil.which", side_effect=which_without_rna_tools),
            patch("petra.cli.launch_workflow", return_value=0) as launch,
        ):
            result = CliRunner().invoke(app, self.run_arguments())

        self.assertEqual(result.exit_code, 0, result.output)
        launch.assert_called_once()


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
