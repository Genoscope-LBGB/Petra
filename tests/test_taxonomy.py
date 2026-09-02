from __future__ import annotations

import subprocess
import unittest
from unittest.mock import patch

from petra.taxonomy import HostTaxonError, LINEAGE_RANKS, resolve_host_taxon

_LINEAGE_ARGUMENTS = ["lineage", "--no-lineage", "--show-name", "--show-rank"]


class ResolveHostTaxonTest(unittest.TestCase):
    def test_resolves_a_taxid_with_a_single_taxonkit_call(self) -> None:
        with patch(
            "petra.taxonomy.subprocess.run",
            side_effect=[_process("9606\tHomo sapiens\tspecies\n")],
        ) as run:
            result = resolve_host_taxon("9606")

        self.assertEqual(result, "Homo sapiens")
        self.assertEqual(run.call_count, 1)
        self.assertEqual(
            run.call_args_list[0].args[0], ["taxonkit", *_LINEAGE_ARGUMENTS]
        )

    def test_resolves_a_name_with_two_taxonkit_calls(self) -> None:
        with patch(
            "petra.taxonomy.subprocess.run",
            side_effect=[
                _process("human\t9606\n"),
                _process("9606\tHomo sapiens\tspecies\n"),
            ],
        ) as run:
            result = resolve_host_taxon("human", taxonkit_executable="custom-taxonkit")

        self.assertEqual(result, "Homo sapiens")
        self.assertEqual(run.call_count, 2)
        self.assertEqual(
            run.call_args_list[0].args[0],
            ["custom-taxonkit", "name2taxid"],
        )
        self.assertEqual(run.call_args_list[0].kwargs["input"], "human\n")
        self.assertEqual(
            run.call_args_list[1].args[0],
            ["custom-taxonkit", *_LINEAGE_ARGUMENTS],
        )

    def test_accepts_a_domain_however_ncbi_spells_its_rank(self) -> None:
        for rank in ("domain", "superkingdom", "acellular root"):
            with self.subTest(rank=rank):
                with patch(
                    "petra.taxonomy.subprocess.run",
                    side_effect=[_process(f"2\tBacteria\t{rank}\n")],
                ):
                    self.assertEqual(resolve_host_taxon("2"), "Bacteria")

    def test_rejects_an_unknown_name(self) -> None:
        with patch(
            "petra.taxonomy.subprocess.run",
            return_value=_process("unknown\t\n"),
        ):
            with self.assertRaisesRegex(
                HostTaxonError,
                "Pass an NCBI TaxID or a scientific name",
            ):
                resolve_host_taxon("unknown")

    def test_rejects_an_ambiguous_name_and_lists_candidate_taxids(self) -> None:
        with patch(
            "petra.taxonomy.subprocess.run",
            return_value=_process("Drosophila\t7215\nDrosophila\t32281\n"),
        ):
            with self.assertRaises(HostTaxonError) as caught:
                resolve_host_taxon("Drosophila")

        self.assertIn("7215", str(caught.exception))
        self.assertIn("32281", str(caught.exception))
        self.assertIn("Pass one of these TaxIDs", str(caught.exception))

    def test_rejects_a_taxon_at_a_rank_petra_does_not_report(self) -> None:
        with patch(
            "petra.taxonomy.subprocess.run",
            side_effect=[_process("1912919\tNasutitermitinae\tsubtribe\n")],
        ):
            with self.assertRaises(HostTaxonError) as caught:
                resolve_host_taxon("1912919")

        self.assertIn("is a subtribe", str(caught.exception))
        for rank in LINEAGE_RANKS:
            self.assertIn(rank, str(caught.exception))

    def test_reports_a_taxonkit_failure_with_stderr(self) -> None:
        with patch(
            "petra.taxonomy.subprocess.run",
            return_value=_process("", stderr="taxonomy database missing", returncode=1),
        ):
            with self.assertRaises(HostTaxonError) as caught:
                resolve_host_taxon("9606")

        self.assertIn("taxonkit lineage failed with exit code 1", str(caught.exception))
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
