import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "src" / "petra" / "workflow" / "scripts"
)
sys.path.insert(0, str(SCRIPTS_DIRECTORY))

from filter_overrepresented_targets import filter_alignment  # noqa: E402


class FilterAlignmentTest(unittest.TestCase):
    def test_removes_every_alignment_to_targets_above_the_maximum(self) -> None:
        retained_at_limit = [
            f"query-{index}\tat-limit\t1e-20\t100\t1\n" for index in range(50)
        ]
        removed_above_limit = [
            f"query-{index}\tover-limit\t1e-20\t100\t2\n" for index in range(51)
        ]
        retained_single = "query-single\tsingle\t1e-20\t100\t3\n"

        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            source = directory / "diamond.tsv"
            destination = directory / "filtered.tsv"
            source.write_text(
                "".join(retained_at_limit)
                + "".join(removed_above_limit)
                + retained_single,
                encoding="utf-8",
            )

            filter_alignment(source, destination, maximum=50)

            self.assertEqual(
                destination.read_text(encoding="utf-8"),
                "".join(retained_at_limit) + retained_single,
            )


if __name__ == "__main__":
    unittest.main()
