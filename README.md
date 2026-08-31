# PETRA - Predicted protein-enabled Taxonomic Rank Assignment

## Installation

```console
pip install .
```

## Usage

### Running the workflow

Supply every long-read file by repeating `--long-reads`. The sequencing technology
defaults to `pacbio` (`map-hifi`), select `ont` to use `map-ont`:

```console
petra run \
  --assembly assembly.fasta \
  --diamond-db database/petra.dmnd \
  --long-reads reads-1.fastq.gz \
  --long-reads reads-2.fastq.gz \
  --technology ont \
  --output run-1
```

The run writes `predictions.summary.png`, a summary graph to eyeball the result:
GC against depth of coverage with point area scaled by sequence length, both
cumulative-length marginals, the five longest sequences, and the share of
sequences and bases that received an assignment.

Customize it with `--summary-font-size`, `--summary-disable-panel` (repeat to
drop `header`, `gc-histogram`, `coverage-histogram`, `taxonomy-legend`,
`longest-sequences`, `length-legend` or `assignment-status`) and
`--summary-colors`. Sequences are colored by superkingdom unless
`--summary-colors` supplies a two-column TSV of taxonomic rank and color:

```text
Actinopteri	#1F78B4
Bacteria	green
```

A sequence then takes the color of the deepest assigned rank listed in the
file, assigned sequences matching no rank are drawn as `other`, and unassigned
sequences keep the unassigned color. The same graph can be rendered from an
existing table without rerunning the workflow:

```console
python plot_summary_graph.py -i run-1/predictions.final.tsv -o summary.png
```

### Comparing runs

`petra-compare`puts the blobplot of several
finished runs side by side to compare how
assignments landed. Repeat `--directory` once per PETRA output directory and
`--name` once per organism, in the same order:

```console
petra-compare \
  --directory run-1 --name "Chelon auratus" \
  --directory run-2 --name "Pachymatisma johnstonia" \
  --directory run-3 --name "Jania squamata" \
  --colors taxon_colors.tsv \
  --font-size 12 \
  --output comparison.png
```

### Building the local DIAMOND database

Build a taxonomy-aware DIAMOND database from local inputs:

```console
petra build-db \
  --fasta proteins.fasta \
  --taxon-map prot.accession2taxid.FULL.gz \
  --output database
```

The command requires `diamond` on `PATH` and creates `database/petra.dmnd`. If
an input is omitted, it downloads NCBI NR or the full NCBI protein taxon map
into the output directory, logs the source and destination, and removes the
download after the build. It fails rather than replacing an existing download
or `petra.dmnd`.
