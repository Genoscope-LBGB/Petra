from __future__ import annotations

import subprocess
import unittest
from pathlib import Path
from unittest.mock import patch

from petra.taxonomy import HostTaxonError, resolve_host_phylum



class ResolveHostPhylumTest(unittest.TestCase):
    taxonomy = Path("database/taxonomy")

    def test_resolves_a_taxid(self) -> None:
        with patch(
            "petra.taxonomy.subprocess.run",
            return_value=_process("9606\tChordata\n"),
        ):
            self.assertEqual(
                resolve_host_phylum("9606", taxonomy=self.taxonomy), "Chordata"
            )

    def test_resolves_a_name(self) -> None:
        with patch(
            "petra.taxonomy.subprocess.run",
            side_effect=[
                _process("human\t9606\n"),
                _process("9606\tChordata\n"),
            ],
        ):
            self.assertEqual(
                resolve_host_phylum("human", taxonomy=self.taxonomy), "Chordata"
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

    def test_returns_no_phylum_for_an_unranked_taxon(self) -> None:
        with patch(
            "petra.taxonomy.subprocess.run",
            return_value=_process("3451528\tN/A\n"),
        ):
            self.assertIsNone(resolve_host_phylum("3451528", taxonomy=self.taxonomy))

    def test_rejects_an_unknown_taxid(self) -> None:
        with (
            patch(
                "petra.taxonomy.subprocess.run",
                return_value=_process("999999999\tUNKNOWN_TAXID\n"),
            ),
            self.assertRaises(HostTaxonError),
        ):
            resolve_host_phylum("999999999", taxonomy=self.taxonomy)

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
