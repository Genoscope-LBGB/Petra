"""Render the PETRA summary graph from the final taxonomy table."""

import argparse
import csv
import math
import traceback
from collections.abc import Iterable, Sequence as SequenceABC
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.axes import Axes  # noqa: E402
from matplotlib.colors import is_color_like, to_rgba  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from matplotlib.legend import Legend  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.ticker import FuncFormatter, NullFormatter  # noqa: E402

from assign_taxonomy import TAXONOMIC_RANKS  # noqa: E402

_MISSING_TAXON = "N/A"
_LOW_CONFIDENCE_REASON = "low_confidence"
_NO_CONSENSUS_REASON = "no_consensus"
_WHOKARYOTE_REASON = "whokaryote"
_METADATA_COLUMNS = (
    "sequence_name",
    "sequence_length",
    "gc_percent",
    "median_coverage",
    "number_of_aligned_clusters",
)
_REQUIRED_COLUMNS = (*_METADATA_COLUMNS, *TAXONOMIC_RANKS, "unassigned_reason")

PANELS = (
    "header",
    "gc-histogram",
    "coverage-histogram",
    "taxonomy-legend",
    "longest-sequences",
    "length-legend",
    "assignment-status",
)
_LEGEND_PANELS = (
    "taxonomy-legend",
    "longest-sequences",
    "length-legend",
    "assignment-status",
)

_SUPERKINGDOM_COLORS = {
    "Bacteria": "#00CD00",
    "Archaea": "#0000CD",
    "Viruses": "#CD00CD",
    "Eukaryota": "#B22222",
}
_UNASSIGNED_GROUP = "N/A"
_UNASSIGNED_COLOR = "#FFFFFF"
_OTHER_GROUP = "other"
_OTHER_COLOR = "#808080"
_EDGE_COLOR = "#3A3A3A"
_PANEL_BACKGROUND = "#EBEBEB"
_GRID_COLOR = "#FFFFFF"
_MUTED_COLOR = "#555555"
_FILL_ALPHA = 0.55

_ASSIGNED_STATUS = "assigned"
_STATUS_COLORS = {
    _ASSIGNED_STATUS: "#2166AC",
    "low confidence": "#8AB8D8",
    "no consensus": "#FDB863",
    "whokaryote veto": "#B79ACC",
    "no protein evidence": "#D9D9D9",
    "unassigned (other)": "#9E9E9E",
}

_LONGEST_MARKERS = ("o", "^", "s", "P", "X")
_MINIMUM_MARKER_AREA = 18.0
_MAXIMUM_MARKER_AREA = 700.0
_LONGEST_MARKER_AREA = 80.0
_HISTOGRAM_BINS = 33
_LEGEND_LEFT = 0.03
_LEGEND_TOP_WITH_HEADER = 0.825
_LEGEND_TOP = 0.935
_LEGEND_BOTTOM = 0.105
_MINIMUM_LEGEND_GAP = 0.005
_LEGEND_GAP_LINES = 3.0
_MINIMUM_LEGEND_SCALE = 0.4
_LEGEND_FIT_ATTEMPTS = 4
_STATUS_BAR_HEIGHT = 0.05
_STATUS_BAR_PADDING = 0.012
_STATUS_TICK_WIDTH = 0.026
_STATUS_TICK_FONT = 8.5
_STATUS_TITLE_HEIGHT = 0.022
_STATUS_LABEL_THRESHOLD = 12.0


@dataclass(frozen=True)
class AssemblySequence:
    """One row of the final taxonomy table."""

    name: str
    length: int
    gc_percent: float
    median_coverage: float
    aligned_clusters: int
    lineage: tuple[str, ...]
    reason: str

    @property
    def assigned(self) -> bool:
        return any(taxon != _MISSING_TAXON for taxon in self.lineage)


@dataclass(frozen=True)
class Fonts:
    """Every text size in the figure, derived from one base size."""

    base: float

    @property
    def title(self) -> float:
        return self.base * 2.1

    @property
    def header(self) -> float:
        return self.base * 0.95

    @property
    def legend_title(self) -> float:
        return self.base * 1.2

    @property
    def legend(self) -> float:
        return self.base

    @property
    def longest_legend(self) -> float:
        return self.base * 0.9

    @property
    def axis_label(self) -> float:
        return self.base * 1.1

    @property
    def marginal_label(self) -> float:
        return self.base

    @property
    def tick(self) -> float:
        return self.base * 0.9

    @property
    def status_legend(self) -> float:
        return self.base * 0.75

    @property
    def status_tick(self) -> float:
        return self.base * 0.85

    @property
    def status_value(self) -> float:
        return self.base * 0.8


def _read_sequences(source: Path) -> list[AssemblySequence]:
    sequences: list[AssemblySequence] = []
    with source.open(encoding="utf-8", newline="") as input_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        missing = [
            column
            for column in _REQUIRED_COLUMNS
            if column not in (reader.fieldnames or ())
        ]
        if missing:
            raise ValueError(
                f"{source} is missing required columns: {', '.join(missing)}"
            )
        for line_number, row in enumerate(reader, start=2):
            try:
                sequences.append(
                    AssemblySequence(
                        name=row["sequence_name"],
                        length=int(row["sequence_length"]),
                        gc_percent=float(row["gc_percent"]),
                        median_coverage=float(row["median_coverage"]),
                        aligned_clusters=int(row["number_of_aligned_clusters"]),
                        lineage=tuple(row[rank] for rank in TAXONOMIC_RANKS),
                        reason=row["unassigned_reason"],
                    )
                )
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"{source} holds an invalid value on line {line_number}: {error}"
                ) from error
    if not sequences:
        raise ValueError(f"{source} does not hold any assembly sequence.")
    return sequences


def _read_color_map(source: Path) -> dict[str, str]:
    """Read a two-column taxon/color TSV, tolerating one header line."""
    colors: dict[str, str] = {}
    with source.open(encoding="utf-8", newline="") as input_handle:
        for line_number, row in enumerate(
            csv.reader(input_handle, delimiter="\t"), start=1
        ):
            if not row or not row[0].strip() or row[0].lstrip().startswith("#"):
                continue
            if len(row) < 2:
                raise ValueError(
                    f"{source} needs a taxon and a color column on line {line_number}."
                )
            taxon, color = row[0].strip(), row[1].strip()
            if not is_color_like(color):
                if line_number == 1 and not colors:
                    continue
                raise ValueError(
                    f"{source} holds an unknown color on line {line_number}: {color!r}"
                )
            if taxon == _MISSING_TAXON:
                raise ValueError(
                    f"{source} cannot color {_MISSING_TAXON!r}: unassigned sequences "
                    "always keep the unassigned color."
                )
            if taxon in colors:
                raise ValueError(
                    f"{source} holds a duplicate taxon on line {line_number}: {taxon!r}"
                )
            colors[taxon] = color
    if not colors:
        raise ValueError(f"{source} does not hold any taxon color.")
    return colors


def _group_of(sequence: AssemblySequence, color_map: dict[str, str] | None) -> str:
    """Return the color group: superkingdom by default, deepest mapped taxon otherwise."""
    if not sequence.assigned:
        return _UNASSIGNED_GROUP
    if color_map is None:
        domain = sequence.lineage[0]
        return domain if domain != _MISSING_TAXON else _OTHER_GROUP
    for taxon in reversed(sequence.lineage):
        if taxon != _MISSING_TAXON and taxon in color_map:
            return taxon
    return _OTHER_GROUP


def _palette(
    groups: SequenceABC[str], color_map: dict[str, str] | None
) -> tuple[list[str], dict[str, str]]:
    used = set(groups)
    if color_map is None:
        ordered = [group for group in _SUPERKINGDOM_COLORS if group in used]
        palette = {group: _SUPERKINGDOM_COLORS[group] for group in ordered}
    else:
        ordered = [group for group in color_map if group in used]
        palette = {group: color_map[group] for group in ordered}
    for group, color in ((_OTHER_GROUP, _OTHER_COLOR), (_UNASSIGNED_GROUP, _UNASSIGNED_COLOR)):
        if group in used:
            ordered.append(group)
            palette[group] = color
    return ordered, palette


def _status_of(sequence: AssemblySequence) -> str:
    if sequence.assigned:
        return _ASSIGNED_STATUS
    if sequence.reason == _LOW_CONFIDENCE_REASON:
        return "low confidence"
    if sequence.reason == _NO_CONSENSUS_REASON:
        return "no consensus"
    if sequence.reason == _WHOKARYOTE_REASON:
        return "whokaryote veto"
    if sequence.aligned_clusters == 0:
        return "no protein evidence"
    return "unassigned (other)"


def _deepest_assignment(sequence: AssemblySequence) -> str:
    for rank, taxon in zip(reversed(TAXONOMIC_RANKS), reversed(sequence.lineage)):
        if taxon != _MISSING_TAXON:
            return f"{taxon} ({rank})"
    return "unassigned"


def _bp_unit(value: float) -> tuple[str, float]:
    for name, divisor in (("Gbp", 1e9), ("Mbp", 1e6), ("kbp", 1e3)):
        if value >= divisor:
            return name, divisor
    return "bp", 1.0


def _human_bp(value: float) -> str:
    name, divisor = _bp_unit(value)
    if divisor == 1.0:
        return f"{value:.0f} bp"
    return f"{f'{value / divisor:.2f}'.rstrip('0').rstrip('.')} {name}"


def _round_significant(value: float) -> float:
    if value <= 0:
        return 0.0
    magnitude = 10 ** math.floor(math.log10(value))
    return max(magnitude, round(value / magnitude) * magnitude)


def _n50(lengths_longest_first: SequenceABC[int], total: int) -> int:
    cumulative = 0
    for length in lengths_longest_first:
        cumulative += length
        if cumulative >= total / 2:
            return length
    return lengths_longest_first[-1]


def _weighted_median(pairs: Iterable[tuple[float, float]]) -> float:
    ordered = sorted(pairs)
    half = sum(weight for _, weight in ordered) / 2
    cumulative = 0.0
    for value, weight in ordered:
        cumulative += weight
        if cumulative >= half:
            return value
    return ordered[-1][0]


def _linear_bins(lowest: float, highest: float, count: int) -> list[float]:
    if highest <= lowest:
        lowest, highest = lowest - 1.0, highest + 1.0
    step = (highest - lowest) / count
    return [lowest + index * step for index in range(count + 1)]


def _logarithmic_bins(lowest: float, highest: float, count: int) -> list[float]:
    lowest = max(lowest, 1e-6)
    highest = max(highest, lowest * 10)
    start, step = math.log10(lowest), (math.log10(highest) - math.log10(lowest)) / count
    return [10 ** (start + index * step) for index in range(count + 1)]


def _area_scale(lengths: SequenceABC[int]):
    roots = [math.sqrt(length) for length in lengths]
    lowest, span = min(roots), max(roots) - min(roots)

    def scale(value: float) -> float:
        if span <= 0:
            return _MAXIMUM_MARKER_AREA
        ratio = min(max((math.sqrt(value) - lowest) / span, 0.0), 1.0)
        return _MINIMUM_MARKER_AREA + ratio * (
            _MAXIMUM_MARKER_AREA - _MINIMUM_MARKER_AREA
        )

    return scale


def _length_legend_values(lengths: SequenceABC[int]) -> list[float]:
    lowest, highest = max(min(lengths), 1), max(lengths)
    if highest <= lowest:
        return [_round_significant(highest)]
    values = [
        _round_significant(lowest * (highest / lowest) ** (index / 3))
        for index in range(4)
    ]
    return sorted(dict.fromkeys(values))


def _style_panel(axis: Axes, fonts: Fonts, *, grid_axis: str = "both") -> None:
    axis.set_facecolor(_PANEL_BACKGROUND)
    for spine in axis.spines.values():
        spine.set_visible(False)
    axis.grid(True, color=_GRID_COLOR, linewidth=0.8, axis=grid_axis, zorder=0)
    axis.grid(
        True,
        color=_GRID_COLOR,
        linewidth=0.35,
        which="minor",
        axis=grid_axis,
        zorder=0,
    )
    axis.tick_params(length=0, labelsize=fonts.tick)


def _draw_blobplot(
    axis: Axes,
    sequences: SequenceABC[AssemblySequence],
    groups: SequenceABC[str],
    palette: dict[str, str],
    scale,
    fonts: Fonts,
    *,
    show_longest: bool,
) -> None:
    plotted = [
        index
        for index in sorted(
            range(len(sequences)), key=lambda index: -sequences[index].length
        )
        if sequences[index].median_coverage > 0
    ]
    # Largest first so small sequences are never buried under a long one.
    axis.scatter(
        [sequences[index].gc_percent for index in plotted],
        [sequences[index].median_coverage for index in plotted],
        s=[scale(sequences[index].length) for index in plotted],
        marker="o",
        facecolor=[to_rgba(palette[groups[index]], _FILL_ALPHA) for index in plotted],
        edgecolor=_EDGE_COLOR,
        linewidth=0.6,
        zorder=3,
    )
    if show_longest:
        longest = [
            index
            for index in range(min(len(_LONGEST_MARKERS), len(sequences)))
            if sequences[index].median_coverage > 0
        ]
        for rank in reversed(longest):
            axis.scatter(
                [sequences[rank].gc_percent],
                [sequences[rank].median_coverage],
                s=_LONGEST_MARKER_AREA,
                marker=_LONGEST_MARKERS[rank],
                facecolor="black",
                edgecolor="white",
                linewidth=1.1,
                zorder=6 + (len(_LONGEST_MARKERS) - rank),
            )
    coverages = [sequences[index].median_coverage for index in plotted]
    axis.set_yscale("log")
    axis.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    axis.yaxis.set_minor_formatter(NullFormatter())
    axis.set_xlabel("GC (%)", fontsize=fonts.axis_label)
    axis.set_ylabel("Depth of coverage (X)", fontsize=fonts.axis_label)
    _style_panel(axis, fonts)
    if coverages:
        axis.set_ylim(min(coverages) * 0.7, max(coverages) * 1.7)


def _stacked_lengths(
    sequences: SequenceABC[AssemblySequence],
    groups: SequenceABC[str],
    ordered_groups: SequenceABC[str],
    key,
) -> tuple[list[list[float]], list[list[float]]]:
    values: list[list[float]] = []
    weights: list[list[float]] = []
    for group in ordered_groups:
        selected = [
            index
            for index in range(len(sequences))
            if groups[index] == group and sequences[index].median_coverage > 0
        ]
        values.append([key(sequences[index]) for index in selected])
        weights.append([float(sequences[index].length) for index in selected])
    return values, weights


def _draw_gc_histogram(
    axis: Axes,
    sequences: SequenceABC[AssemblySequence],
    groups: SequenceABC[str],
    ordered_groups: SequenceABC[str],
    palette: dict[str, str],
    bins: list[float],
    unit: tuple[str, float],
    fonts: Fonts,
) -> None:
    values, weights = _stacked_lengths(
        sequences, groups, ordered_groups, lambda sequence: sequence.gc_percent
    )
    axis.hist(
        values,
        bins=bins,
        weights=weights,
        stacked=True,
        color=[palette[group] for group in ordered_groups],
        edgecolor="none",
        zorder=3,
    )
    axis.set_ylabel(f"Cumulative\nlength ({unit[0]})", fontsize=fonts.marginal_label)
    axis.tick_params(labelbottom=False)
    axis.yaxis.set_major_formatter(
        FuncFormatter(lambda value, _: f"{value / unit[1]:g}")
    )
    _style_panel(axis, fonts)


def _draw_coverage_histogram(
    axis: Axes,
    sequences: SequenceABC[AssemblySequence],
    groups: SequenceABC[str],
    ordered_groups: SequenceABC[str],
    palette: dict[str, str],
    bins: list[float],
    unit: tuple[str, float],
    fonts: Fonts,
) -> None:
    values, weights = _stacked_lengths(
        sequences, groups, ordered_groups, lambda sequence: sequence.median_coverage
    )
    axis.hist(
        values,
        bins=bins,
        weights=weights,
        stacked=True,
        orientation="horizontal",
        color=[palette[group] for group in ordered_groups],
        edgecolor="none",
        zorder=3,
    )
    axis.set_xlabel(f"Cumulative\nlength ({unit[0]})", fontsize=fonts.marginal_label)
    axis.tick_params(labelleft=False)
    axis.xaxis.set_major_formatter(
        FuncFormatter(lambda value, _: f"{value / unit[1]:g}")
    )
    axis.yaxis.set_minor_formatter(NullFormatter())
    _style_panel(axis, fonts)


def _taxonomy_handles(
    sequences: SequenceABC[AssemblySequence],
    groups: SequenceABC[str],
    ordered_groups: SequenceABC[str],
    palette: dict[str, str],
) -> list[Line2D]:
    handles: list[Line2D] = []
    for group in ordered_groups:
        selected = [
            sequences[index] for index in range(len(sequences)) if groups[index] == group
        ]
        bases = sum(sequence.length for sequence in selected)
        handles.append(
            Line2D(
                [],
                [],
                marker="o",
                linestyle="none",
                markersize=8,
                markerfacecolor=to_rgba(palette[group], _FILL_ALPHA),
                markeredgecolor=_EDGE_COLOR,
                markeredgewidth=0.6,
                label=f"{group}  ({len(selected)} seq, {_human_bp(bases)})",
            )
        )
    return handles


def _longest_handles(
    sequences: SequenceABC[AssemblySequence],
    groups: SequenceABC[str],
    palette: dict[str, str],
) -> list[Line2D]:
    return [
        Line2D(
            [],
            [],
            marker=_LONGEST_MARKERS[rank],
            linestyle="none",
            markersize=8,
            markerfacecolor=palette[groups[rank]],
            markeredgecolor=_EDGE_COLOR,
            markeredgewidth=0.8,
            label=(
                f"{sequences[rank].name} - {_human_bp(sequences[rank].length)} - "
                f"{sequences[rank].gc_percent:.1f}% GC - "
                f"{sequences[rank].median_coverage:g}X\n"
                f"      {_deepest_assignment(sequences[rank])}"
            ),
        )
        for rank in range(min(len(_LONGEST_MARKERS), len(sequences)))
    ]


def _length_handles(values: SequenceABC[float], scale) -> list[Line2D]:
    return [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markersize=math.sqrt(scale(value)),
            markerfacecolor="white",
            markeredgecolor=_EDGE_COLOR,
            label=_human_bp(value),
        )
        for value in values
    ]


def _draw_status_bars(
    axis: Axes,
    statuses: SequenceABC[str],
    ordered_statuses: SequenceABC[str],
    lengths: SequenceABC[int],
    fonts: Fonts,
) -> None:
    counts = [statuses.count(status) for status in ordered_statuses]
    bases = [
        sum(lengths[index] for index in range(len(statuses)) if statuses[index] == status)
        for status in ordered_statuses
    ]
    rows = (
        (bases, float(sum(lengths)), "bases"),
        ([float(count) for count in counts], float(len(statuses)), "seq."),
    )
    for row, (values, total, _) in enumerate(rows):
        left = 0.0
        for status, value in zip(ordered_statuses, values):
            percent = 100 * value / total if total else 0.0
            axis.barh(
                row,
                percent,
                left=left,
                color=_STATUS_COLORS[status],
                edgecolor="white",
                height=0.62,
            )
            if percent > _STATUS_LABEL_THRESHOLD:
                axis.text(
                    left + percent / 2,
                    row,
                    f"{percent:.0f}%",
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontsize=fonts.status_value,
                    color="white" if status == _ASSIGNED_STATUS else "#333333",
                    weight="bold",
                )
            left += percent
    axis.set_yticks([0, 1], [label for _, _, label in rows], fontsize=fonts.status_tick)
    axis.set_xlim(0, 100)
    axis.set_xticks([0, 50, 100])
    axis.tick_params(length=0, labelsize=fonts.status_value)
    for spine in axis.spines.values():
        spine.set_visible(False)


@dataclass(frozen=True)
class _LegendBlock:
    """One titled legend block of the right-hand column."""

    title: str
    handles: list
    font_size: float
    labels: list[str] | None = None
    extra_height: float = 0.0
    keywords: dict = field(default_factory=dict)


def _fit_text_width(figure: Figure, artist, available: float) -> None:
    """Shrink a text artist until it fits the width it was given."""
    inverted = figure.transFigure.inverted()
    for _ in range(12):
        figure.canvas.draw()
        width = inverted.transform_bbox(artist.get_window_extent()).width
        if width <= available or artist.get_fontsize() <= 4.0:
            return
        artist.set_fontsize(artist.get_fontsize() * max(available / width, 0.8))


def _create_legend_blocks(
    axis_legend: Axes, blocks: SequenceABC[_LegendBlock], fonts: Fonts, scale: float
) -> list[Legend]:
    created: list[Legend] = []
    for block in blocks:
        keywords = dict(block.keywords)
        if block.labels is not None:
            keywords["labels"] = block.labels
        legend = axis_legend.legend(
            handles=block.handles,
            title=block.title,
            loc="upper left",
            bbox_to_anchor=(_LEGEND_LEFT, 0.5),
            frameon=False,
            title_fontproperties={"size": fonts.legend_title * scale},
            alignment="left",
            fontsize=block.font_size * scale,
            **keywords,
        )
        axis_legend.add_artist(legend)
        created.append(legend)
    return created


def _draw_legend_column(
    figure: Figure,
    axis_legend: Axes,
    blocks: SequenceABC[_LegendBlock],
    fonts: Fonts,
    top: float,
) -> float:
    """Lay the column out with one shared gap, shrinking fonts only if it overflows.

    The last block may own the assignment-status bars through ``extra_height``; it stays
    anchored to the bottom so its legend never drifts away from its bars.
    """
    position = axis_legend.get_position()
    available_height = top - _LEGEND_BOTTOM
    available_width = position.width * (1.0 - _LEGEND_LEFT)
    scale = 1.0
    for _ in range(_LEGEND_FIT_ATTEMPTS):
        created = _create_legend_blocks(axis_legend, blocks, fonts, scale)
        figure.canvas.draw()
        inverted = figure.transFigure.inverted()
        boxes = [inverted.transform_bbox(item.get_window_extent()) for item in created]
        heights = [
            box.height + block.extra_height for box, block in zip(boxes, blocks)
        ]
        reserved = _MINIMUM_LEGEND_GAP * max(len(created) - 1, 0)
        overflow = max(
            sum(heights) / max(available_height - reserved, 1e-6),
            max(box.width for box in boxes) / available_width,
            1.0,
        )
        if overflow <= 1.0 or scale <= _MINIMUM_LEGEND_SCALE:
            break
        for item in created:
            item.remove()
        scale = max(scale / overflow, _MINIMUM_LEGEND_SCALE)

    bottom_anchored = blocks[-1].extra_height > 0
    bottom = _LEGEND_BOTTOM + (heights[-1] if bottom_anchored else 0.0)
    movable = created[:-1] if bottom_anchored else list(created)
    movable_heights = heights[:-1] if bottom_anchored else list(heights)
    gap_count = len(movable) if bottom_anchored else max(len(movable) - 1, 1)
    even_gap = (top - bottom - sum(movable_heights)) / gap_count
    # Cap the gap so a sparse column packs at the top instead of drifting apart.
    cap = _LEGEND_GAP_LINES * fonts.base * scale / 72.0 / figure.get_figheight()
    gap = min(max(even_gap, _MINIMUM_LEGEND_GAP), cap)

    def anchor(legend: Legend, block_top: float) -> None:
        legend.set_bbox_to_anchor(
            (_LEGEND_LEFT, (block_top - position.y0) / position.height),
            transform=axis_legend.transAxes,
        )

    offset = top
    for item, height in zip(movable, movable_heights):
        anchor(item, offset)
        offset -= height + gap
    if bottom_anchored:
        anchor(created[-1], bottom)
    return scale


def plot_summary_graph(
    taxonomy: Path,
    destination: Path,
    *,
    colors: Path | None = None,
    disabled_panels: Iterable[str] = (),
    font_size: float = 10.0,
    width: float = 16.5,
    height: float = 9.6,
    dpi: int = 130,
    title: str = "PETRA summary graph",
) -> None:
    disabled = tuple(disabled_panels)
    unknown = [panel for panel in disabled if panel not in PANELS]
    if unknown:
        raise ValueError(
            f"Unknown panel: {', '.join(sorted(unknown))}. "
            f"Available panels: {', '.join(PANELS)}."
        )
    if font_size <= 0:
        raise ValueError("Font size must be greater than zero.")
    enabled = {panel for panel in PANELS if panel not in disabled}
    fonts = Fonts(base=font_size)

    color_map = _read_color_map(colors) if colors is not None else None
    sequences = sorted(_read_sequences(taxonomy), key=lambda item: -item.length)
    groups = [_group_of(sequence, color_map) for sequence in sequences]
    ordered_groups, palette = _palette(groups, color_map)
    statuses = [_status_of(sequence) for sequence in sequences]
    ordered_statuses = [status for status in _STATUS_COLORS if status in set(statuses)]

    lengths = [sequence.length for sequence in sequences]
    total = sum(lengths)
    unit = _bp_unit(float(total))
    scale = _area_scale(lengths)
    covered = [sequence for sequence in sequences if sequence.median_coverage > 0]
    hidden = len(sequences) - len(covered)
    if not covered:
        raise ValueError(f"{taxonomy} does not hold any sequence with coverage above 0.")
    median_coverage = _weighted_median(
        (sequence.median_coverage, float(sequence.length)) for sequence in covered
    )

    show_gc = "gc-histogram" in enabled
    show_coverage = "coverage-histogram" in enabled
    legend_panels = [panel for panel in _LEGEND_PANELS if panel in enabled]
    show_column = bool(legend_panels) or "header" in enabled

    figure = plt.figure(figsize=(width, height), dpi=dpi)
    figure.patch.set_facecolor("white")
    width_ratios = [3.0] + ([0.95] if show_coverage else []) + ([1.45] if show_column else [])
    height_ratios = ([1.0] if show_gc else []) + [3.3]
    grid = GridSpec(
        len(height_ratios),
        len(width_ratios),
        figure=figure,
        width_ratios=width_ratios,
        height_ratios=height_ratios,
        hspace=0.06,
        wspace=0.06,
        left=0.075,
        right=0.985,
        top=0.955,
        bottom=0.10,
    )
    main_row = 1 if show_gc else 0
    axis_main = figure.add_subplot(grid[main_row, 0])
    _draw_blobplot(
        axis_main,
        sequences,
        groups,
        palette,
        scale,
        fonts,
        show_longest="longest-sequences" in enabled,
    )

    if show_gc:
        axis_gc = figure.add_subplot(grid[0, 0], sharex=axis_main)
        _draw_gc_histogram(
            axis_gc,
            sequences,
            groups,
            ordered_groups,
            palette,
            _linear_bins(
                math.floor(min(sequence.gc_percent for sequence in covered)) - 1,
                math.ceil(max(sequence.gc_percent for sequence in covered)) + 1,
                _HISTOGRAM_BINS,
            ),
            unit,
            fonts,
        )
    if show_coverage:
        axis_coverage = figure.add_subplot(grid[main_row, 1], sharey=axis_main)
        _draw_coverage_histogram(
            axis_coverage,
            sequences,
            groups,
            ordered_groups,
            palette,
            _logarithmic_bins(
                min(sequence.median_coverage for sequence in covered) * 0.7,
                max(sequence.median_coverage for sequence in covered) * 1.7,
                _HISTOGRAM_BINS,
            ),
            unit,
            fonts,
        )

    if show_column:
        axis_legend = figure.add_subplot(grid[:, len(width_ratios) - 1])
        axis_legend.axis("off")
        position = axis_legend.get_position()
        align_x = position.x0 + _LEGEND_LEFT * position.width

        if "header" in enabled:
            header_lines = (
                (0.962, title, fonts.title, "normal"),
                (
                    0.925,
                    f"{len(sequences)} sequences | {_human_bp(total)} | "
                    f"N50 {_human_bp(_n50(lengths, total))} | "
                    f"median {median_coverage:g}X",
                    fonts.header,
                    "normal",
                ),
                (
                    0.895,
                    f"({hidden} 0X sequences have been hidden)",
                    fonts.legend,
                    "italic",
                ),
            )
            for offset, text, size, style in header_lines:
                artist = axis_legend.text(
                    0.5,
                    offset,
                    text,
                    horizontalalignment="center",
                    fontsize=size,
                    style=style,
                    color="black" if size == fonts.title else _MUTED_COLOR,
                    transform=axis_legend.transAxes,
                )
                _fit_text_width(figure, artist, position.width)

        blocks: list[_LegendBlock] = []
        if "taxonomy-legend" in enabled:
            blocks.append(
                _LegendBlock(
                    title="taxon" if color_map is not None else "superkingdom",
                    handles=_taxonomy_handles(
                        sequences, groups, ordered_groups, palette
                    ),
                    font_size=fonts.legend,
                )
            )
        if "longest-sequences" in enabled:
            blocks.append(
                _LegendBlock(
                    title="Longest sequences",
                    handles=_longest_handles(sequences, groups, palette),
                    font_size=fonts.longest_legend,
                    keywords={"labelspacing": 0.8},
                )
            )
        if "length-legend" in enabled:
            blocks.append(
                _LegendBlock(
                    title="Length",
                    handles=_length_handles(_length_legend_values(lengths), scale),
                    font_size=fonts.legend,
                    keywords={"labelspacing": 1.3, "handletextpad": 1.5},
                )
            )
        if "assignment-status" in enabled:
            blocks.append(
                _LegendBlock(
                    title="Assignment status",
                    handles=[
                        Rectangle((0, 0), 1, 1, color=_STATUS_COLORS[status])
                        for status in ordered_statuses
                    ],
                    font_size=fonts.status_legend,
                    labels=[
                        f"{status} ({statuses.count(status)})"
                        for status in ordered_statuses
                    ],
                    extra_height=(
                        _STATUS_BAR_PADDING + _STATUS_BAR_HEIGHT + _STATUS_TITLE_HEIGHT
                    ),
                    keywords={
                        "ncol": 2,
                        "handlelength": 0.9,
                        "columnspacing": 0.9,
                        "handletextpad": 0.5,
                    },
                )
            )

        if blocks:
            _draw_legend_column(
                figure,
                axis_legend,
                blocks,
                fonts,
                _LEGEND_TOP_WITH_HEADER if "header" in enabled else _LEGEND_TOP,
            )
        if "assignment-status" in enabled:
            tick_width = _STATUS_TICK_WIDTH * fonts.status_tick / _STATUS_TICK_FONT
            axis_status = figure.add_axes(
                [
                    align_x + tick_width,
                    _LEGEND_BOTTOM,
                    position.x1 - (align_x + tick_width),
                    _STATUS_BAR_HEIGHT,
                ]
            )
            _draw_status_bars(axis_status, statuses, ordered_statuses, lengths, fonts)

    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, facecolor="white")
    plt.close(figure)
    print(
        f"Rendered {destination} from {len(sequences)} assembly sequences "
        f"({_human_bp(total)}) with panels: {', '.join(sorted(enabled))}."
    )


def _parse_arguments(argv: SequenceABC[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--taxonomy",
        type=Path,
        required=True,
        help="Final PETRA taxonomy table (predictions.final.tsv).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Image to write; the extension selects the format.",
    )
    parser.add_argument(
        "--colors",
        type=Path,
        help=(
            "Two-column TSV of taxon and color. A sequence is colored with the "
            "deepest assigned rank found in the file; assigned sequences without a "
            "match are colored as 'other'. Superkingdom colors are used when omitted."
        ),
    )
    parser.add_argument(
        "--disable-panel",
        action="append",
        choices=PANELS,
        default=[],
        metavar="PANEL",
        help=f"Panel to leave out; repeatable. Choices: {', '.join(PANELS)}.",
    )
    parser.add_argument(
        "--font-size", type=float, default=10.0, help="Base font size in points."
    )
    parser.add_argument("--width", type=float, default=16.5, help="Figure width in inches.")
    parser.add_argument(
        "--height", type=float, default=9.6, help="Figure height in inches."
    )
    parser.add_argument("--dpi", type=int, default=130, help="Figure resolution.")
    parser.add_argument(
        "--title", default="PETRA summary graph", help="Title shown in the header."
    )
    return parser.parse_args(argv)


def main(argv: SequenceABC[str] | None = None) -> None:
    arguments = _parse_arguments(argv)
    plot_summary_graph(
        arguments.taxonomy,
        arguments.output,
        colors=arguments.colors,
        disabled_panels=arguments.disable_panel,
        font_size=arguments.font_size,
        width=arguments.width,
        height=arguments.height,
        dpi=arguments.dpi,
        title=arguments.title,
    )


if "snakemake" in globals():
    workflow = globals()["snakemake"]
    log_path = Path(workflow.log[0])
    with log_path.open("w", encoding="utf-8") as log_handle:
        with redirect_stdout(log_handle), redirect_stderr(log_handle):
            try:
                color_inputs = list(workflow.input.colors)
                plot_summary_graph(
                    Path(workflow.input.taxonomy),
                    Path(workflow.output.graph),
                    colors=Path(color_inputs[0]) if color_inputs else None,
                    disabled_panels=workflow.params.disabled_panels,
                    font_size=float(workflow.params.font_size),
                )
            except Exception:
                traceback.print_exc()
                raise
elif __name__ == "__main__":
    main()
