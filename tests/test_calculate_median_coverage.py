import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "src" / "petra" / "workflow" / "scripts"
)
sys.path.insert(0, str(SCRIPTS_DIRECTORY))

from calculate_median_coverage import calculate_median_coverage  # noqa: E402


class CalculateMedianCoverageTest(unittest.TestCase):
    def test_weights_depth_by_interval_length(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            source = directory / "coverage.bedgraph"
            destination = directory / "coverage.tsv"
            source.write_text(
                "contig-a\t0\t3\t0\n"
                "contig-a\t3\t4\t10\n"
                "contig-b\t0\t1\t1\n"
                "contig-b\t1\t2\t2\n"
                "contig-c\t0\t5\t7\n",
                encoding="utf-8",
            )

            calculate_median_coverage(source, destination)

            self.assertEqual(
                destination.read_text(encoding="utf-8"),
                "sequence_name\tmedian_coverage\n"
                "contig-a\t0\n"
                "contig-b\t1.5\n"
                "contig-c\t7\n",
            )


if __name__ == "__main__":
    unittest.main()
