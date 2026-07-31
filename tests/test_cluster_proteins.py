import csv
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
    def test_reports_member_count_for_each_retained_cluster(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            source = directory / "proteins.faa"
            clustered = directory / "clustered.faa"
            cluster_sizes = directory / "cluster_sizes.tsv"
            source.write_text(
                ">contig_1\nM\n>contig_2\nM\n>contig_3\nM\n",
                encoding="utf-8",
            )

            def run_cdhit(group: ProteinGroup) -> Path:
                group.clustered.write_text(
                    ">contig_1 representative\nM\n>contig_3\nM\n",
                    encoding="utf-8",
                )
                Path(f"{group.clustered}.clstr").write_text(
                    ">Cluster 0\n"
                    "0 1aa, >contig_1 representative... *\n"
                    "1 1aa, >contig_2... at 100.00%\n"
                    ">Cluster 1\n"
                    "0 1aa, >contig_3... *\n",
                    encoding="utf-8",
                )
                return group.clustered

            cluster_proteins(
                source,
                clustered,
                cluster_sizes,
                workers=1,
                max_clusters=10,
                runner=run_cdhit,
            )

            self.assertEqual(
                clustered.read_text(encoding="utf-8"),
                ">contig_1 representative\nM\n>contig_3\nM\n",
            )
            with cluster_sizes.open(encoding="utf-8", newline="") as input_handle:
                rows = list(csv.reader(input_handle, delimiter="\t"))
            self.assertEqual(
                rows,
                [
                    ["cluster_representative", "number_of_genes_in_cluster"],
                    ["contig_1", "2"],
                    ["contig_3", "1"],
                ],
            )


if __name__ == "__main__":
    unittest.main()
