import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "src" / "petra" / "workflow" / "scripts"
)
sys.path.insert(0, str(SCRIPTS_DIRECTORY))

from assign_taxonomy import TAXONOMIC_RANKS  # noqa: E402
from plot_summary_comparison import (  # noqa: E402
    _UNASSIGNED_LABEL,
    _legend_handles,
    _load_run,
    _panel_letter,
    _taxonomy_file,
    main,
    plot_summary_comparison,
)
from plot_summary_graph import _palette  # noqa: E402

_COLUMNS = (
    "sequence_name",
    "sequence_length",
    "gc_percent",
    "median_coverage",
    "number_of_aligned_unigenes",
    *TAXONOMIC_RANKS,
    "unassigned_reason",
)


def _row(
    name: str,
    length: int,
    gc_percent: float,
    coverage: float,
    lineage: dict[str, str],
    reason: str = "-",
) -> str:
    values = {
        "sequence_name": name,
        "sequence_length": str(length),
        "gc_percent": f"{gc_percent}",
        "median_coverage": f"{coverage:g}",
        "number_of_aligned_unigenes": "12",
        "unassigned_reason": reason,
        **{rank: lineage.get(rank, "N/A") for rank in TAXONOMIC_RANKS},
    }
    return "\t".join(values[column] for column in _COLUMNS)


def _run_directory(directory: Path, name: str, rows: tuple[str, ...]) -> Path:
    destination = directory / name
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "predictions.final.tsv").write_text(
        "\n".join(("\t".join(_COLUMNS), *rows, "")), encoding="utf-8"
    )
    return destination


_FISH = (
    _row("contig-1", 30_000_000, 42.0, 46, {"domain": "Eukaryota", "phylum": "Chordata"}),
    _row("contig-2", 68_000, 49.8, 8, {"domain": "Eukaryota"}),
    _row("contig-3", 21_000, 46.0, 2, {}, "insufficient_support"),
)
_SPONGE = (
    _row("contig-1", 900_000, 60.1, 25, {"domain": "Bacteria"}),
    _row("contig-2", 400_000, 38.0, 90, {"domain": "Archaea"}),
    _row("contig-3", 15_000, 44.5, 3, {}, "no_consensus"),
)


class PanelLetterTest(unittest.TestCase):
    def test_counts_past_the_alphabet(self) -> None:
        self.assertEqual(
            [_panel_letter(index) for index in (0, 1, 25, 26, 27, 51, 52)],
            ["a", "b", "z", "aa", "ab", "az", "ba"],
        )


class TaxonomyFileTest(unittest.TestCase):
    def test_resolves_a_run_directory_to_its_table(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = _run_directory(Path(temporary_directory), "fish", _FISH)

            self.assertEqual(
                _taxonomy_file(directory), directory / "predictions.final.tsv"
            )

    def test_accepts_the_table_itself(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            table = (
                _run_directory(Path(temporary_directory), "fish", _FISH)
                / "predictions.final.tsv"
            )

            self.assertEqual(_taxonomy_file(table), table)

    def test_rejects_a_directory_without_a_table(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            with self.assertRaisesRegex(ValueError, "predictions.final.tsv"):
                _taxonomy_file(Path(temporary_directory))


class SharedPaletteTest(unittest.TestCase):
    def test_gives_one_color_per_taxon_across_every_run(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            runs = [
                _load_run("Fish", _run_directory(directory, "fish", _FISH), None),
                _load_run("Sponge", _run_directory(directory, "sponge", _SPONGE), None),
            ]
            ordered, palette = _palette(
                [group for run in runs for group in run.groups], None
            )

            self.assertEqual(ordered, ["Bacteria", "Archaea", "Eukaryota", "N/A"])
            self.assertEqual(len(set(palette.values())), len(ordered))

    def test_labels_unassigned_sequences_in_the_legend(self) -> None:
        ordered, palette = _palette(["Bacteria", "N/A"], None)
        labels = [handle.get_label() for handle in _legend_handles(ordered, palette)]

        self.assertEqual(labels, ["Bacteria", _UNASSIGNED_LABEL])


class LoadRunTest(unittest.TestCase):
    def test_rejects_a_table_without_a_covered_sequence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = _run_directory(
                Path(temporary_directory),
                "flat",
                (_row("contig-1", 1_000, 42.0, 0, {"domain": "Bacteria"}),),
            )

            with self.assertRaisesRegex(ValueError, "coverage above 0"):
                _load_run("Flat", directory, None)

    def test_rejects_an_empty_table(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = _run_directory(Path(temporary_directory), "empty", ())

            with self.assertRaisesRegex(ValueError, "any assembly sequence"):
                _load_run("Empty", directory, None)


class PlotSummaryComparisonTest(unittest.TestCase):
    def test_renders_one_panel_per_run(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            image = directory / "comparison.png"
            plot_summary_comparison(
                [
                    ("Chelon auratus", _run_directory(directory, "fish", _FISH)),
                    ("Pachymatisma", _run_directory(directory, "sponge", _SPONGE)),
                ],
                image,
            )

            self.assertGreater(image.stat().st_size, 10_000)

    def test_renders_with_a_shared_color_file_and_a_large_font(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            colors = directory / "colors.tsv"
            colors.write_text("Chordata\t#FF0000\nBacteria\tgreen\n", encoding="utf-8")
            image = directory / "comparison.png"
            plot_summary_comparison(
                [
                    ("Chelon auratus", _run_directory(directory, "fish", _FISH)),
                    ("Pachymatisma", _run_directory(directory, "sponge", _SPONGE)),
                ],
                image,
                colors=colors,
                font_size=18.0,
            )

            self.assertGreater(image.stat().st_size, 10_000)

    def test_rejects_a_font_size_of_zero(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            with self.assertRaisesRegex(ValueError, "Font size"):
                plot_summary_comparison(
                    [("Fish", _run_directory(directory, "fish", _FISH))],
                    directory / "comparison.png",
                    font_size=0.0,
                )

    def test_rejects_an_empty_run_list(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            with self.assertRaisesRegex(ValueError, "At least one run"):
                plot_summary_comparison([], Path(temporary_directory) / "empty.png")


class MainTest(unittest.TestCase):
    def test_pairs_every_directory_with_a_name(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            image = directory / "comparison.png"
            main(
                [
                    "-d",
                    str(_run_directory(directory, "fish", _FISH)),
                    "-n",
                    "Chelon auratus",
                    "-d",
                    str(_run_directory(directory, "sponge", _SPONGE)),
                    "-n",
                    "Pachymatisma johnstonia",
                    "-o",
                    str(image),
                ]
            )

            self.assertTrue(image.is_file())

    def test_rejects_more_directories_than_names(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            with self.assertRaisesRegex(ValueError, "one name per directory"):
                main(
                    [
                        "-d",
                        str(_run_directory(directory, "fish", _FISH)),
                        "-d",
                        str(_run_directory(directory, "sponge", _SPONGE)),
                        "-n",
                        "Only one",
                        "-o",
                        str(directory / "comparison.png"),
                    ]
                )


if __name__ == "__main__":
    unittest.main()
