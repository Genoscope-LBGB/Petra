from __future__ import annotations

import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import yaml

from petra.database import TAXONOMY_DIRECTORY_NAME, TAXONOMY_FILENAMES
from petra.launcher import workflow_directory


class WorkflowDagTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.assembly = self.directory / "assembly.fasta"
        self.reads = self.directory / "reads.fastq"
        self.diamond_database = self.directory / "petra.dmnd"
        self.covariance_model = self.directory / "petra.rrna.cm"
        self.clan_info = self.directory / "petra.rrna.clanin"
        self.ssu_database = self.directory / "petra.ssu.fasta"
        self.lsu_database = self.directory / "petra.lsu.fasta"
        self.taxonomy = self.directory / TAXONOMY_DIRECTORY_NAME
        self.taxonomy.mkdir()
        for filename in TAXONOMY_FILENAMES:
            (self.taxonomy / filename).touch()
        for path in (
            self.assembly,
            self.reads,
            self.diamond_database,
            self.covariance_model,
            self.clan_info,
            self.ssu_database,
            self.lsu_database,
        ):
            path.touch()

    def test_optional_analyses_only_enter_the_dag_when_enabled(self) -> None:
        protein_only = self._dry_run(rrna=False, trna=False)
        self._assert_job_count(protein_only, "split_assembly_by_domain", 1)
        self._assert_job_count(protein_only, "extract_organelles", 1)
        self._assert_job_absent(protein_only, "cmscan")
        self._assert_job_absent(protein_only, "summarize_rrna_hits")
        self._assert_job_absent(protein_only, "trnascan")
        self.assertIn("predictions.organelles.tsv", protein_only)
        self.assertIn("predictions.gene_taxonomy.tsv", protein_only)
        self.assertIn("whokaryote.noneuk.txt", protein_only)
        for domain in (
            "eukaryota",
            "bacteria",
            "viruses",
            "archaea",
            "unassigned",
        ):
            self.assertIn(
                f"decontaminated/predictions.{domain}.fasta",
                protein_only,
            )

        rrna = self._dry_run(rrna=True, trna=False)
        self._assert_job_count(rrna, "split_assembly_by_domain", 1)
        self._assert_job_count(rrna, "cmscan", 16)
        self._assert_job_count(rrna, "summarize_rrna_hits", 1)
        self._assert_job_absent(rrna, "trnascan")
        self.assertIn("rna/predictions.rrna.gff", rrna)
        self.assertIn("rna/predictions.rrna.tsv", rrna)

        trna = self._dry_run(rrna=False, trna=True)
        self._assert_job_count(trna, "split_assembly_by_domain", 1)
        self._assert_job_absent(trna, "cmscan")
        self._assert_job_absent(trna, "summarize_rrna_hits")
        self._assert_job_count(trna, "trnascan", 1)
        self.assertIn("rna/predictions.trna.tsv", trna)
        self.assertIn("rna/predictions.trna.stats", trna)

        rna_only = self._dry_run(
            rrna=True,
            trna=True,
            skip_taxonomic_assignment=True,
        )
        self._assert_job_absent(rna_only, "split_assembly_by_domain")
        self._assert_job_absent(rna_only, "diamond")
        self._assert_job_absent(rna_only, "extract_organelles")
        self._assert_job_absent(rna_only, "plot_summary_graph")
        self.assertNotIn("predictions.gene_taxonomy.tsv", rna_only)
        self._assert_job_count(rna_only, "cmscan", 16)
        self._assert_job_count(rna_only, "summarize_rrna_hits", 1)
        self._assert_job_count(rna_only, "trnascan", 1)

    def _dry_run(
        self,
        *,
        rrna: bool,
        trna: bool,
        skip_taxonomic_assignment: bool = False,
    ) -> str:
        config = {
            "assembly": str(self.assembly),
            "diamond_db": str(self.diamond_database),
            "long_reads": [str(self.reads)],
            "sequencing_technology": "pacbio",
            "no_bam_cleaning": False,
            "max_low_complexity": 100,
            "max_clusters_per_sequence": 1000,
            "max_target_occurrences": 50,
            "alpha": 0.5,
            "min_aligned_clusters": 3,
            "cluster_support_length": 50_000,
            "clusters_per_decade": 3,
            "summary_font_size": 10.0,
            "summary_colors": None,
            "summary_disabled_panels": [],
            "rrna": rrna,
            "trna": trna,
            "skip_taxonomic_assignment": skip_taxonomic_assignment,
            "rrna_chunks": 16,
            "rrna_min_alignment_length": 200,
            "host_phylum": None,
            "covariance_model": str(self.covariance_model),
            "clan_info": str(self.clan_info),
            "ssu_database": str(self.ssu_database),
            "lsu_database": str(self.lsu_database),
            "taxonomy": str(self.taxonomy),
        }
        config_path = self.directory / (
            f"config-{rrna}-{trna}-{skip_taxonomic_assignment}.yaml"
        )
        config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
        process = subprocess.run(
            [
                str(Path(sys.executable).parent / "snakemake"),
                "--snakefile",
                str(workflow_directory() / "Snakefile"),
                "--configfile",
                str(config_path),
                "--directory",
                str(self.directory),
                "--cores",
                "1",
                "--dry-run",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        output = process.stdout + process.stderr
        self.assertEqual(process.returncode, 0, output)
        return output

    def _assert_job_count(self, output: str, rule: str, count: int) -> None:
        self.assertRegex(output, rf"(?m)^{re.escape(rule)}\s+{count}\s*$")

    def _assert_job_absent(self, output: str, rule: str) -> None:
        self.assertNotRegex(output, rf"(?m)^{re.escape(rule)}\s+\d+\s*$")


if __name__ == "__main__":
    unittest.main()
