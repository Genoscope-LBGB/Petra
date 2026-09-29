import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIRECTORY = (
    Path(__file__).resolve().parents[1] / "src" / "petra" / "workflow" / "scripts"
)
sys.path.insert(0, str(SCRIPTS_DIRECTORY))

from assign_taxonomy import TAXONOMIC_RANKS  # noqa: E402
from plot_summary_graph import (  # noqa: E402
    PANELS,
    THEMES,
    _compact_bp,
    _group_of,
    _human_bp,
    _palette,
    _read_color_map,
    _read_sequences,
    _row_major,
    plot_summary_graph,
    theme_colors,
)

_COLUMNS = (
    "sequence_name",
    "sequence_length",
    "gc_percent",
    "median_coverage",
    "number_of_aligned_unigenes",
    *TAXONOMIC_RANKS,
    "unassigned_reason",
)
_ACTINOPTERI = {
    "domain": "Eukaryota",
    "kingdom": "Metazoa",
    "phylum": "Chordata",
    "class": "Actinopteri",
}


def _row(
    name: str,
    length: int,
    gc_percent: float,
    coverage: float,
    unigenes: int,
    lineage: dict[str, str],
    reason: str = "-",
) -> str:
    values = {
        "sequence_name": name,
        "sequence_length": str(length),
        "gc_percent": f"{gc_percent}",
        "median_coverage": f"{coverage:g}",
        "number_of_aligned_unigenes": str(unigenes),
        "unassigned_reason": reason,
        **{rank: lineage.get(rank, "N/A") for rank in TAXONOMIC_RANKS},
    }
    return "\t".join(values[column] for column in _COLUMNS)


def _taxonomy(directory: Path) -> Path:
    destination = directory / "predictions.final.tsv"
    destination.write_text(
        "\n".join(
            (
                "\t".join(_COLUMNS),
                _row("contig-1", 30_000_000, 42.0, 46, 259, _ACTINOPTERI),
                _row("contig-2", 68_000, 49.8, 8, 12, {"domain": "Eukaryota", "kingdom": "Metazoa"}),
                _row("contig-3", 21_000, 46.0, 2, 17, {}, "no_consensus"),
                _row("contig-4", 19_500, 54.7, 59, 2, {}, "insufficient_support"),
                _row("contig-5", 13_100, 50.7, 127, 0, {}),
                "",
            )
        ),
        encoding="utf-8",
    )
    return destination


class ColorGroupTest(unittest.TestCase):
    def test_colors_a_sequence_with_its_deepest_mapped_rank(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            sequences = _read_sequences(_taxonomy(directory))
            color_map = {"Chordata": "#111111", "Actinopteri": "#222222"}

            self.assertEqual(_group_of(sequences[0], color_map), "Actinopteri")

    def test_falls_back_to_other_when_no_rank_is_mapped(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            sequences = _read_sequences(_taxonomy(directory))

            self.assertEqual(_group_of(sequences[1], {"Actinopteri": "#222222"}), "other")

    def test_keeps_unassigned_sequences_unassigned(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            sequences = _read_sequences(_taxonomy(directory))

            self.assertEqual(_group_of(sequences[2], {"Chordata": "#111111"}), "N/A")

    def test_atlasea_colors_eukaryotes_without_a_listed_phylum_as_eukaryota(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            sequences = _read_sequences(_taxonomy(directory))
            color_map = _read_color_map(theme_colors("atlasea"))

            self.assertEqual(_group_of(sequences[0], color_map), "Chordata")
            self.assertEqual(_group_of(sequences[1], color_map), "Eukaryota")


class ThemeTest(unittest.TestCase):
    def test_every_packaged_theme_is_a_valid_color_map(self) -> None:
        self.assertIn("superkingdom", THEMES)
        self.assertIn("atlasea", THEMES)
        for theme in THEMES:
            with self.subTest(theme=theme):
                _read_color_map(theme_colors(theme))

    def test_rejects_an_unknown_theme(self) -> None:
        with self.assertRaisesRegex(ValueError, "Unknown theme: rainbow"):
            theme_colors("rainbow")


class ColorMapTest(unittest.TestCase):
    def _write(self, directory: Path, content: str) -> Path:
        destination = directory / "colors.tsv"
        destination.write_text(content, encoding="utf-8")
        return destination

    def test_skips_a_header_line_and_keeps_file_order(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            source = self._write(
                Path(temporary_directory),
                "taxon\tcolor\nActinopteri\t#1F78B4\nBacteria\tgreen\n",
            )

            self.assertEqual(
                _read_color_map(source),
                {"Actinopteri": "#1F78B4", "Bacteria": "green"},
            )

    def test_rejects_an_unknown_color(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            source = self._write(
                Path(temporary_directory), "Actinopteri\t#1F78B4\nBacteria\tnotacolor\n"
            )

            with self.assertRaisesRegex(ValueError, "unknown color"):
                _read_color_map(source)

    def test_rejects_coloring_unassigned_sequences(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            source = self._write(Path(temporary_directory), "N/A\tred\n")

            with self.assertRaisesRegex(ValueError, "unassigned"):
                _read_color_map(source)

    def test_rejects_a_duplicate_taxon(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            source = self._write(
                Path(temporary_directory), "Actinopteri\tred\nActinopteri\tblue\n"
            )

            with self.assertRaisesRegex(ValueError, "duplicate taxon"):
                _read_color_map(source)


class HumanUnitTest(unittest.TestCase):
    def test_scales_to_the_largest_fitting_unit(self) -> None:
        self.assertEqual(_human_bp(900), "900 bp")
        self.assertEqual(_human_bp(11_107), "11.11 kbp")
        self.assertEqual(_human_bp(1_000_000), "1 Mbp")
        self.assertEqual(_human_bp(696_502_509), "696.5 Mbp")
        self.assertEqual(_human_bp(2_400_000_000), "2.4 Gbp")

    def test_compacts_to_three_significant_digits(self) -> None:
        self.assertEqual(_compact_bp(900), "900 bp")
        self.assertEqual(_compact_bp(410_440), "410 kbp")
        self.assertEqual(_compact_bp(53_530_000), "53.5 Mbp")
        self.assertEqual(_compact_bp(999_600), "1 Mbp")


class LegendOrderTest(unittest.TestCase):
    def test_orders_taxa_by_total_length_with_other_and_unassigned_last(self) -> None:
        ordered, _ = _palette(
            ["Bacteria", "N/A", "Eukaryota", "other", "Eukaryota"],
            [300, 9_000, 200, 5_000, 200],
            {"Bacteria": "green", "Eukaryota": "red"},
        )

        self.assertEqual(ordered, ["Eukaryota", "Bacteria", "other", "N/A"])

    def test_fills_legend_columns_so_they_read_row_by_row(self) -> None:
        self.assertEqual(_row_major(["a", "b", "c", "d", "e"], 2), ["a", "c", "e", "b", "d"])


class RenderTest(unittest.TestCase):
    def test_renders_every_panel_by_default(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            destination = directory / "predictions.summary.png"

            plot_summary_graph(_taxonomy(directory), destination)

            self.assertGreater(destination.stat().st_size, 10_000)

    def test_renders_without_any_optional_panel(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            destination = directory / "blobplot.png"

            plot_summary_graph(
                _taxonomy(directory),
                destination,
                disabled_panels=PANELS,
                font_size=6.0,
            )

            self.assertGreater(destination.stat().st_size, 10_000)

    def test_renders_with_a_color_map(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            colors = directory / "colors.tsv"
            colors.write_text("Actinopteri\t#1F78B4\n", encoding="utf-8")
            destination = directory / "predictions.summary.png"

            plot_summary_graph(_taxonomy(directory), destination, colors=colors)

            self.assertGreater(destination.stat().st_size, 10_000)

    def test_rejects_an_unknown_panel(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)

            with self.assertRaisesRegex(ValueError, "Unknown panel: histogram"):
                plot_summary_graph(
                    _taxonomy(directory),
                    directory / "graph.png",
                    disabled_panels=["histogram"],
                )

    def test_rejects_a_table_without_the_required_columns(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            source = directory / "broken.tsv"
            source.write_text("sequence_name\tsequence_length\ncontig-1\t10\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "missing required columns"):
                plot_summary_graph(source, directory / "graph.png")


if __name__ == "__main__":
    unittest.main()
