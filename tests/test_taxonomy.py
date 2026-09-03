from __future__ import annotations

import subprocess
import unittest
from pathlib import Path
from unittest.mock import patch

from petra.taxonomy import HostTaxonError, resolve_host_phylum

_PHYLUM_ARGUMENTS = [
    "reformat2",
    "--taxid-field",
    "1",
    "--format",
    "{phylum}",
    "--miss-rank-repl",
    "N/A",
    "--miss-taxid-repl",
    "N/A",
]


class ResolveHostPhylumTest(unittest.TestCase):
    taxonomy = Path("database/taxonomy")

    def test_resolves_a_taxid_with_a_single_taxonkit_call(self) -> None:
        with patch(
            "petra.taxonomy.subprocess.run",
            side_effect=[_process("9606\tChordata\n")],
        ) as run:
            result = resolve_host_phylum("9606", taxonomy=self.taxonomy)

        self.assertEqual(result, "Chordata")
        self.assertEqual(run.call_count, 1)
        self.assertEqual(
            run.call_args_list[0].args[0],
            ["taxonkit", "--data-dir", str(self.taxonomy), *_PHYLUM_ARGUMENTS],
        )
        self.assertEqual(run.call_args_list[0].kwargs["input"], "9606\n")

    def test_resolves_a_name_with_two_taxonkit_calls(self) -> None:
        with patch(
            "petra.taxonomy.subprocess.run",
            side_effect=[
                _process("human\t9606\n"),
                _process("9606\tChordata\n"),
            ],
        ) as run:
            result = resolve_host_phylum(
                "human",
                taxonomy=self.taxonomy,
                taxonkit_executable="custom-taxonkit",
            )

        self.assertEqual(result, "Chordata")
        self.assertEqual(run.call_count, 2)
        self.assertEqual(
            run.call_args_list[0].args[0],
            [
                "custom-taxonkit",
                "--data-dir",
                str(self.taxonomy),
                "name2taxid",
            ],
        )
        self.assertEqual(run.call_args_list[0].kwargs["input"], "human\n")
        self.assertEqual(
            run.call_args_list[1].args[0],
            [
                "custom-taxonkit",
                "--data-dir",
                str(self.taxonomy),
                *_PHYLUM_ARGUMENTS,
            ],
        )

    def test_accepts_a_taxid_at_phylum_rank(self) -> None:
        with patch(
            "petra.taxonomy.subprocess.run",
            return_value=_process("7711\tChordata\n"),
        ):
            self.assertEqual(
                resolve_host_phylum("7711", taxonomy=self.taxonomy),
                "Chordata",
            )

    def test_rejects_a_taxon_without_a_phylum(self) -> None:
        with (
            patch(
                "petra.taxonomy.subprocess.run",
                return_value=_process("2\tN/A\n"),
            ),
            self.assertRaisesRegex(
                HostTaxonError,
                "TaxID 2 has no phylum in the NCBI taxonomy",
            ),
        ):
            resolve_host_phylum("2", taxonomy=self.taxonomy)

    def test_rejects_an_unknown_name(self) -> None:
        with (
            patch(
                "petra.taxonomy.subprocess.run",
                return_value=_process("unknown\t\n"),
            ),
            self.assertRaisesRegex(
                HostTaxonError,
                "Pass an NCBI TaxID or a scientific name",
            ),
        ):
            resolve_host_phylum("unknown", taxonomy=self.taxonomy)

    def test_rejects_an_ambiguous_name_and_lists_candidate_taxids(self) -> None:
        with (
            patch(
                "petra.taxonomy.subprocess.run",
                return_value=_process("Drosophila\t7215\nDrosophila\t32281\n"),
            ),
            self.assertRaises(HostTaxonError) as caught,
        ):
            resolve_host_phylum("Drosophila", taxonomy=self.taxonomy)

        self.assertIn("7215", str(caught.exception))
        self.assertIn("32281", str(caught.exception))
        self.assertIn("Pass one of these TaxIDs", str(caught.exception))

    def test_reports_a_taxonkit_failure_with_stderr(self) -> None:
        with (
            patch(
                "petra.taxonomy.subprocess.run",
                return_value=_process(
                    "", stderr="taxonomy database missing", returncode=1
                ),
            ),
            self.assertRaises(HostTaxonError) as caught,
        ):
            resolve_host_phylum("9606", taxonomy=self.taxonomy)

        self.assertIn(
            "taxonkit reformat2 failed with exit code 1", str(caught.exception)
        )
        self.assertIn("taxonomy database missing", str(caught.exception))


def _process(
    stdout: str,
    *,
    stderr: str = "",
    returncode: int = 0,
) -> subprocess.CompletedProcess[str]:
    return subprocess.CompletedProcess(
        args=["taxonkit"],
        returncode=returncode,
        stdout=stdout,
        stderr=stderr,
    )


if __name__ == "__main__":
    unittest.main()
