# PETRA - Predicted protein-enabled Taxonomic Rank Assignment

## Installation

```console
pip install .
```

## Usage

### Building the local databases

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

Add `--build-rna-dbs` to also build the databases used for rRNA detection, next
to `petra.dmnd` in the same output directory:

- `petra.rrna.cm`, the `cmpress`ed concatenation of the 5S, 5.8S, SSU and LSU
  Rfam families, for `cmscan`.
- `petra.rrna.clanin`, built alongside `petra.rrna.cm` and filtered to the
  families present in that CM, for `cmscan` clan competition.
- `petra.ssu.fasta`, a BLAST database of SILVA SSURef NR99 and PR2, carrying
  NCBI taxon identifiers.
- `petra.lsu.fasta`, a BLAST database of SILVA LSURef NR99 and EUKARYOME,
  carrying the lineage in each sequence header.

This part requires `makeblastdb`, `cmpress` and `taxonkit` on `PATH`, plus an
NCBI taxonomy dump that `taxonkit` can read. SILVA sequences of unidentified
organisms are discarded, and PR2 sequences already present in SILVA are dropped
so the SSU database holds each accession once. 
Sources are downloaded into the output directory and removed after the build.

Because the DIAMOND build takes several hours, the two halves can be run
separately. `--skip-diamond-db` leaves an existing `petra.dmnd` untouched:

```console
petra build-db --output database --build-rna-dbs --skip-diamond-db
```

### Running the workflow

For the full workflow, `--database` takes the directory built by
`petra build-db`, not an individual database file. Supply every long-read file
by repeating `--long-reads`. The sequencing technology defaults to `pacbio`
(`map-hifi`); select `ont` to use `map-ont`:

```console
petra run \
  --assembly assembly.fasta \
  --database database \
  --long-reads reads-1.fastq.gz \
  --long-reads reads-2.fastq.gz \
  --technology ont \
  --output run-1
```

#### Optional rRNA and tRNA analyses

Add `--rrna` to detect rRNA genes with `cmscan` and classify SSU and LSU genes
with `blastn`. It requires Infernal's `cmscan` and BLAST's `blastn` on `PATH`
and produces these files under the output directory's `rna/` subdirectory:

- `rna/predictions.rrna.gff`, the full rRNA gene inventory, including 5S and
  5.8S.
- `rna/predictions.rrna.tsv`, the best BLAST hit per SSU or LSU gene, with a
  pipe-joined lineage and a `non_host` column.

`--host-taxon` accepts an NCBI scientific name or TaxID. PETRA resolves the
taxon to its NCBI phylum and matches that phylum as a whole lineage field.
Without a host, `non_host` is `-`. A hit whose lineage does not resolve to the
host phylum is reported as non-host. Use `--rrna-min-alignment-length` to change
the minimum reported BLAST alignment length from its default of 200.

Add `--trna` to annotate tRNAs. It requires `tRNAscan-SE` on `PATH` and writes
`rna/predictions.trna.tsv` without pseudogenes, undetermined tRNAs or suppressor
tRNAs, plus `rna/predictions.trna.stats`.

Use `--skip-taxonomic-assignment` with `--rrna`, `--trna`, or both to run only
the selected RNA analyses. Long reads are not required in this mode. The
database is still required for rRNA classification, but a tRNA-only run does
not require one:

```console
petra run \
  --assembly assembly.fasta \
  --database database \
  --rrna \
  --trna \
  --skip-taxonomic-assignment \
  --output run-1
```

A full run also writes `predictions.organelles.tsv`, listing each sequence that
Tiara classified as an organelle and its second-stage type (`mitochondrion`,
`plastid`, or `unknown`). As with the other Whokaryote results, only sequences
of at least 5,000 bp are considered.

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

### Filtering the prediction table

`petra-filter` keeps the sequences of `predictions.final.tsv` that satisfy an
expression, given as one quoted argument. A condition is a column name of the
table, an operator among `==`, `!=`, `<`, `<=`, `>`, `>=` (`=` is read as
`==`), and a value:

```console
petra-filter \
  -i run-1/predictions.final.tsv \
  -o kept.tsv \
  "median_coverage<40 and gc_percent>59 and not phylum==Actinopteri"
```

Keywords are case-insensitive, and everything inside the quoted expression, 
parentheses included, reaches the program untouched:

```console
petra-filter -i run-1/predictions.final.tsv \
  "(domain==Bacteria or domain==N/A) and sequence_length>=10000"
```

The filtered table goes to standard output when `--output` is omitted, and the
counts of retained sequences and bases are reported on standard error.
