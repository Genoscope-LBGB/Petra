"""Compose one blobplot per PETRA run so assignments can be compared side by side.

Panels share the taxon palette and one legend; every axis, marker scale and
cumulative-length unit stays per run, because assemblies being compared rarely
share a genome size. The drawing helpers come from ``plot_summary_graph`` so the
two figures cannot drift apart.
"""

import argparse
import sys
from collections.abc import Sequence as SequenceABC
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.legend import Legend  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from matplotlib.text import Text  # noqa: E402
from matplotlib.transforms import Bbox  # noqa: E402

# Snakemake runs these scripts as files, so they import each other flat. The
# installed entry point is imported as a module instead and has to say where
# the siblings live.
_SCRIPTS_DIRECTORY = str(Path(__file__).resolve().parent)
if _SCRIPTS_DIRECTORY not in sys.path:
    sys.path.insert(0, _SCRIPTS_DIRECTORY)

from plot_summary_graph import (  # noqa: E402
    DEFAULT_THEME,
    THEMES,
    AssemblySequence,
    Fonts,
    _EDGE_COLOR,
    _FILL_ALPHA,
    _HISTOGRAM_BINS,
    _UNASSIGNED_GROUP,
    _area_scale,
    _bp_unit,
    _draw_blobplot,
    _draw_coverage_histogram,
    _draw_gc_histogram,
    _group_of,
    _linear_bins,
    _logarithmic_bins,
    _palette,
    _read_color_map,
    _read_sequences,
    _row_major,
    theme_colors,
    to_rgba,
)

_TAXONOMY_FILE = "predictions.final.tsv"
_UNASSIGNED_LABEL = "No assignment"
_NAME_SCALE = 1.7
_LETTER_SCALE = 2.0
_PANEL_WIDTH_RATIOS = (3.2, 1.0)
_PANEL_HEIGHT_RATIOS = (1.0, 3.0)
_PANEL_SPACING = 0.08
_GROUP_SPACING = 0.4
_MARGIN = 0.012
_LEGEND_PADDING = 0.02
_NAME_PADDING = 0.01
_FIT_ATTEMPTS = 3
_LETTERS = "abcdefghijklmnopqrstuvwxyz"


@dataclass(frozen=True)
class ComparedRun:
    """One PETRA run to draw as a panel."""

    name: str
    source: Path
    sequences: tuple[AssemblySequence, ...]
    groups: tuple[str, ...]

    @property
    def covered(self) -> tuple[AssemblySequence, ...]:
        """Sequences the blobplot can place, which needs a positive coverage."""
        return tuple(
            sequence for sequence in self.sequences if sequence.median_coverage > 0
        )


def _panel_letter(index: int) -> str:
    """Return the panel label: a, b, ... z, aa, ab, ..."""
    letter = ""
    while True:
        index, remainder = divmod(index, len(_LETTERS))
        letter = _LETTERS[remainder] + letter
        if index == 0:
            return letter
        index -= 1


def _taxonomy_file(path: Path) -> Path:
    """Resolve a run directory, or an already resolved table, to its table."""
    source = path / _TAXONOMY_FILE if path.is_dir() else path
    if not source.is_file():
        raise ValueError(
            f"{path} does not hold {_TAXONOMY_FILE}: pass a PETRA output directory "
            "or the table itself."
        )
    return source


def _load_run(name: str, path: Path, color_map: dict[str, str]) -> ComparedRun:
    source = _taxonomy_file(path)
    sequences = tuple(_read_sequences(source))
    if not sequences:
        raise ValueError(f"{source} does not hold any assembly sequence.")
    run = ComparedRun(
        name=name,
        source=source,
        sequences=sequences,
        groups=tuple(_group_of(sequence, color_map) for sequence in sequences),
    )
    if not run.covered:
        raise ValueError(f"{source} does not hold any sequence with coverage above 0.")
    return run


def _legend_handles(
    ordered_groups: SequenceABC[str], palette: dict[str, str]
) -> list[Line2D]:
    """One marker per taxon group, without the per-run counts panels cannot share."""
    return [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markersize=9,
            markerfacecolor=to_rgba(palette[group], _FILL_ALPHA),
            markeredgecolor=_EDGE_COLOR,
            markeredgewidth=0.6,
            label=_UNASSIGNED_LABEL if group == _UNASSIGNED_GROUP else group,
        )
        for group in ordered_groups
    ]


def _drawn_bbox(figure: Figure, artist) -> Bbox:
    """Extent of one artist in figure fractions, without forcing a redraw."""
    renderer = figure.canvas.get_renderer()
    measured = (
        artist.get_window_extent(renderer)
        if isinstance(artist, Text)
        else artist.get_tightbbox(renderer)
    )
    return figure.transFigure.inverted().transform_bbox(measured)


def _figure_bbox(figure: Figure, artists) -> Bbox:
    """Union of the drawn extents of artists, in figure fractions."""
    return Bbox.union([_drawn_bbox(figure, artist) for artist in artists])


def _probe_height(figure: Figure, text: str, size: float) -> float:
    """Height a text of that size will need, measured then discarded."""
    probe = figure.text(0.5, 0.5, text, fontsize=size)
    height = _figure_bbox(figure, [probe]).height
    probe.remove()
    return height


def _shrink_to_fit(
    figure: Figure, artist, available: float, *, vertical: bool = False
) -> None:
    """Shrink a text to an extent in one step, as its size scales linearly.

    Iterating here instead would redraw every blobplot per attempt, which costs
    seconds per panel.
    """
    box = _drawn_bbox(figure, artist)
    size = box.height if vertical else box.width
    if size <= available:
        return
    artist.set_fontsize(max(artist.get_fontsize() * available / size, 4.0))


def _fit_axis_labels(figure: Figure, axes) -> None:
    """Keep every axis label inside its own axis, so neighbors cannot collide.

    A large font makes 'Cumulative length' taller than the short marginal it
    belongs to, which is how it ends up drawn over the panel above.
    """
    for axis in axes:
        position = axis.get_position()
        _shrink_to_fit(figure, axis.yaxis.label, position.height, vertical=True)
        _shrink_to_fit(figure, axis.xaxis.label, position.width)


def _draw_legend(figure: Figure, handles: SequenceABC[Line2D], fonts: Fonts) -> Legend:
    """Lay the shared legend on one row, wrapping row by row only if it overflows."""
    available = 1.0 - 2.0 * _MARGIN
    legend = None
    for columns in range(len(handles), 0, -1):
        if legend is not None:
            legend.remove()
        legend = figure.legend(
            handles=_row_major(handles, columns),
            loc="upper center",
            bbox_to_anchor=(0.5, 1.0 - _MARGIN),
            ncol=columns,
            frameon=True,
            fontsize=fonts.legend,
            handletextpad=0.4,
            columnspacing=1.2,
            borderpad=0.6,
        )
        figure.canvas.draw()
        if _figure_bbox(figure, [legend]).width <= available:
            break
    return legend


def plot_summary_comparison(
    runs: SequenceABC[tuple[str, Path]],
    destination: Path,
    *,
    colors: Path | None = None,
    font_size: float = 10.0,
    panel_width: float = 5.5,
    height: float = 5.2,
    dpi: int = 130,
) -> None:
    """Render one comparison figure holding a blobplot per run."""
    if not runs:
        raise ValueError("At least one run is required.")
    for value, label in (
        (font_size, "Font size"),
        (panel_width, "Panel width"),
        (height, "Figure height"),
        (dpi, "Resolution"),
    ):
        if value <= 0:
            raise ValueError(f"{label} must be greater than zero.")

    color_map = _read_color_map(colors if colors is not None else theme_colors(DEFAULT_THEME))
    loaded = [_load_run(name, path, color_map) for name, path in runs]
    ordered_groups, palette = _palette(
        [group for run in loaded for group in run.groups],
        [sequence.length for run in loaded for sequence in run.sequences],
        color_map,
    )

    fonts = Fonts(base=font_size)
    figure = plt.figure(figsize=(panel_width * len(loaded), height), dpi=dpi)
    figure.patch.set_facecolor("white")

    legend = _draw_legend(figure, _legend_handles(ordered_groups, palette), fonts)
    letter_height = _probe_height(
        figure, _panel_letter(len(loaded) - 1), fonts.base * _LETTER_SCALE
    )
    name_height = max(
        _probe_height(figure, run.name, fonts.base * _NAME_SCALE) for run in loaded
    )
    # The letters live between the legend and the panels, so reserve their row.
    top = (
        _figure_bbox(figure, [legend]).y0
        - _LEGEND_PADDING
        - letter_height
        - _NAME_PADDING
    )

    grid = GridSpec(
        1,
        len(loaded),
        figure=figure,
        top=top,
        bottom=name_height + 2.0 * _NAME_PADDING,
        left=_MARGIN,
        right=1.0 - _MARGIN,
        wspace=_GROUP_SPACING,
    )

    groups: list[list] = []
    for index, run in enumerate(loaded):
        cell = grid[0, index].subgridspec(
            2,
            2,
            width_ratios=_PANEL_WIDTH_RATIOS,
            height_ratios=_PANEL_HEIGHT_RATIOS,
            wspace=_PANEL_SPACING,
            hspace=_PANEL_SPACING,
        )
        axis_main = figure.add_subplot(cell[1, 0])
        axis_top = figure.add_subplot(cell[0, 0], sharex=axis_main)
        axis_right = figure.add_subplot(cell[1, 1], sharey=axis_main)
        groups.append([axis_top, axis_main, axis_right])

        covered = run.covered
        total = sum(sequence.length for sequence in run.sequences)
        unit = _bp_unit(float(total))
        scale = _area_scale([sequence.length for sequence in run.sequences])
        _draw_blobplot(
            axis_main,
            run.sequences,
            run.groups,
            palette,
            scale,
            fonts,
            show_longest=False,
        )
        _draw_gc_histogram(
            axis_top,
            run.sequences,
            run.groups,
            ordered_groups,
            palette,
            _linear_bins(
                min(sequence.gc_percent for sequence in covered),
                max(sequence.gc_percent for sequence in covered),
                _HISTOGRAM_BINS,
            ),
            unit,
            fonts,
        )
        _draw_coverage_histogram(
            axis_right,
            run.sequences,
            run.groups,
            ordered_groups,
            palette,
            _logarithmic_bins(
                min(sequence.median_coverage for sequence in covered),
                max(sequence.median_coverage for sequence in covered),
                _HISTOGRAM_BINS,
            ),
            unit,
            fonts,
        )

    _fit_panels(figure, grid, groups, name_height)
    _label_panels(
        figure, grid, groups, [run.name for run in loaded], name_height, fonts
    )

    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, facecolor="white")
    plt.close(figure)
    print(
        f"Rendered {destination} comparing {len(loaded)} runs: "
        f"{', '.join(run.name for run in loaded)}."
    )


def _fit_panels(figure: Figure, grid: GridSpec, groups, name_height: float) -> None:
    """Grow the margins the axis labels and the name row actually need.

    Labels are refitted on every pass: moving a margin resizes the axes they
    have to fit into.
    """
    axes = [axis for group in groups for axis in group]
    needed = name_height + 2.0 * _NAME_PADDING
    for _ in range(_FIT_ATTEMPTS):
        _fit_axis_labels(figure, axes)
        figure.canvas.draw()
        drawn = _figure_bbox(figure, axes)
        left = grid.left + max(_MARGIN - drawn.x0, 0.0)
        right = grid.right - max(drawn.x1 - (1.0 - _MARGIN), 0.0)
        bottom = grid.bottom + max(needed - drawn.y0, 0.0)
        if (left, right, bottom) == (grid.left, grid.right, grid.bottom):
            return
        if right - left <= 0.1 or grid.top - bottom <= 0.1:
            return
        grid.update(left=left, right=right, bottom=bottom)


def _label_panels(
    figure: Figure,
    grid: GridSpec,
    groups,
    names: SequenceABC[str],
    name_height: float,
    fonts: Fonts,
) -> None:
    """Letter every panel and center its name, aligned across the figure.

    Both rows are pinned to the reserved bands rather than to the panel extents,
    which a long axis label overflows.
    """
    figure.canvas.draw()
    boxes = [_figure_bbox(figure, group) for group in groups]
    for index, (box, name) in enumerate(zip(boxes, names)):
        figure.text(
            box.x0,
            grid.top + _NAME_PADDING,
            _panel_letter(index),
            fontsize=fonts.base * _LETTER_SCALE,
            horizontalalignment="left",
            verticalalignment="bottom",
        )
        artist = figure.text(
            (box.x0 + box.x1) / 2.0,
            _NAME_PADDING + name_height / 2.0,
            name,
            fontsize=fonts.base * _NAME_SCALE,
            horizontalalignment="center",
            verticalalignment="center",
        )
        _shrink_to_fit(figure, artist, box.width)


def _parse_arguments(argv: SequenceABC[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-d",
        "--directory",
        action="append",
        required=True,
        type=Path,
        metavar="DIRECTORY",
        help=(
            f"PETRA output directory holding {_TAXONOMY_FILE}, or the table itself; "
            "repeat once per run to compare."
        ),
    )
    parser.add_argument(
        "-n",
        "--name",
        action="append",
        required=True,
        metavar="NAME",
        help="Organism name of each directory, in the same order.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Image to write; the extension selects the format.",
    )
    color_source = parser.add_mutually_exclusive_group()
    color_source.add_argument(
        "--colors",
        type=Path,
        help=(
            "Two-column TSV of taxon and color shared by every panel. A sequence is "
            "colored with the deepest assigned rank found in the file; assigned "
            "sequences without a match are colored as 'other'."
        ),
    )
    color_source.add_argument(
        "--theme",
        choices=THEMES,
        default=DEFAULT_THEME,
        help="Packaged color file to use instead of --colors.",
    )
    parser.add_argument(
        "--font-size",
        type=float,
        default=10.0,
        help="Base font size in points; every other size derives from it.",
    )
    parser.add_argument(
        "--panel-width",
        type=float,
        default=5.5,
        help="Width in inches of one panel; the figure grows with the run count.",
    )
    parser.add_argument(
        "--height", type=float, default=5.2, help="Figure height in inches."
    )
    parser.add_argument("--dpi", type=int, default=130, help="Figure resolution.")
    return parser.parse_args(argv)


def main(argv: SequenceABC[str] | None = None) -> None:
    arguments = _parse_arguments(argv)
    if len(arguments.directory) != len(arguments.name):
        raise ValueError(
            f"{len(arguments.directory)} directories were given for "
            f"{len(arguments.name)} names: pass one name per directory, in order."
        )
    plot_summary_comparison(
        list(zip(arguments.name, arguments.directory)),
        arguments.output,
        colors=arguments.colors or theme_colors(arguments.theme),
        font_size=arguments.font_size,
        panel_width=arguments.panel_width,
        height=arguments.height,
        dpi=arguments.dpi,
    )


def run() -> None:
    """Installed entry point: report a rejected request without a traceback."""
    try:
        main()
    except ValueError as error:
        raise SystemExit(f"Error: {error}")


if __name__ == "__main__":
    run()
