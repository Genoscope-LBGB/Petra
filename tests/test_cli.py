from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from typer.testing import CliRunner

from petra.cli import SequencingTechnology, app


class RunCliTest(unittest.TestCase):
    def test_accepts_repeated_long_reads_and_ont_technology(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            assembly = directory / "assembly.fasta"
            diamond_db = directory / "petra.dmnd"
            first_reads = directory / "first.fastq"
            second_reads = directory / "second.fastq"
            for path in (assembly, diamond_db, first_reads, second_reads):
                path.touch()

            with (
                patch("petra.cli._require_mapping_executables"),
                patch("petra.cli.launch_workflow", return_value=0) as launch,
            ):
                result = CliRunner().invoke(
                    app,
                    [
                        "run",
                        "--assembly",
                        str(assembly),
                        "--diamond-db",
                        str(diamond_db),
                        "--long-reads",
                        str(first_reads),
                        "--long-reads",
                        str(second_reads),
                        "--technology",
                        "ont",
                    ],
                )

            self.assertEqual(result.exit_code, 0, result.output)
            workflow_config = launch.call_args.kwargs["workflow_config"]
            self.assertEqual(workflow_config["long_reads"], [first_reads, second_reads])
            self.assertEqual(
                workflow_config["sequencing_technology"], SequencingTechnology.ONT
            )

    def test_defaults_to_pacbio_technology(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            assembly = directory / "assembly.fasta"
            diamond_db = directory / "petra.dmnd"
            reads = directory / "reads.fastq"
            for path in (assembly, diamond_db, reads):
                path.touch()

            with (
                patch("petra.cli._require_mapping_executables"),
                patch("petra.cli.launch_workflow", return_value=0) as launch,
            ):
                result = CliRunner().invoke(
                    app,
                    [
                        "run",
                        "--assembly",
                        str(assembly),
                        "--diamond-db",
                        str(diamond_db),
                        "--long-reads",
                        str(reads),
                    ],
                )

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(
                launch.call_args.kwargs["workflow_config"]["sequencing_technology"],
                SequencingTechnology.PACBIO,
            )

    def test_rejects_missing_mapping_executables_before_launch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            assembly = directory / "assembly.fasta"
            diamond_db = directory / "petra.dmnd"
            reads = directory / "reads.fastq"
            for path in (assembly, diamond_db, reads):
                path.touch()

            with (
                patch("petra.cli.shutil.which", return_value=None),
                patch("petra.cli.launch_workflow") as launch,
            ):
                result = CliRunner().invoke(
                    app,
                    [
                        "run",
                        "--assembly",
                        str(assembly),
                        "--diamond-db",
                        str(diamond_db),
                        "--long-reads",
                        str(reads),
                    ],
                )

            self.assertEqual(result.exit_code, 1)
            self.assertIn(
                "Required mapping executables not found on PATH: minimap2, samtools.",
                result.output,
            )
            launch.assert_not_called()


if __name__ == "__main__":
    unittest.main()
