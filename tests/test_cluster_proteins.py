import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "src" / "petra" / "workflow" / "scripts"
)
sys.path.insert(0, str(SCRIPTS_DIRECTORY))

from cluster_proteins import ProteinGroup, cluster_proteins  # noqa: E402


class ClusterProteinsTest(unittest.TestCase):
    def test_writes_every_cluster_representative_below_the_limit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            source = directory / "proteins.faa"
            clustered = directory / "clustered.faa"
            source.write_text(
                ">contig_1\nM\n>contig_2\nM\n>contig_3\nM\n",
                encoding="utf-8",
            )

            cluster_proteins(
                source,
                clustered,
                workers=1,
                max_clusters=10,
                runner=self._runner(">contig_1 representative\nM\n>contig_3\nM\n"),
            )

            self.assertEqual(
                clustered.read_text(encoding="utf-8"),
                ">contig_1 representative\nM\n>contig_3\nM\n",
            )

    def test_subsamples_representatives_down_to_the_limit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            source = directory / "proteins.faa"
            clustered = directory / "clustered.faa"
            source.write_text(
                ">contig_1\nM\n>contig_2\nM\n>contig_3\nM\n>contig_4\nM\n",
                encoding="utf-8",
            )

            cluster_proteins(
                source,
                clustered,
                workers=1,
                max_clusters=2,
                runner=self._runner(
                    ">contig_1\nM\n>contig_2\nM\n>contig_3\nM\n>contig_4\nM\n"
                ),
            )

            retained = [
                line[1:]
                for line in clustered.read_text(encoding="utf-8").splitlines()
                if line.startswith(">")
            ]
            self.assertEqual(len(retained), 2)
            self.assertEqual(retained, sorted(retained))
            self.assertLessEqual(
                set(retained),
                {"contig_1", "contig_2", "contig_3", "contig_4"},
            )

    def test_rejects_a_non_positive_cluster_limit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            source = directory / "proteins.faa"
            source.write_text(">contig_1\nM\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "cluster limit must be positive"):
                cluster_proteins(
                    source,
                    directory / "clustered.faa",
                    workers=1,
                    max_clusters=0,
                    runner=self._runner(">contig_1\nM\n"),
                )

    @staticmethod
    def _runner(clustered: str):
        def run_cdhit(group: ProteinGroup) -> Path:
            group.clustered.write_text(clustered, encoding="utf-8")
            return group.clustered

        return run_cdhit


if __name__ == "__main__":
    unittest.main()
