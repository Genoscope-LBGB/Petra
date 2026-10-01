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

The command requires `diamond` on `PATH` and creates `database/petra.dmnd`. It
also downloads the current NCBI taxonomy dump and keeps the TaxonKit files under
`database/taxonomy/`. Pass `--taxdump taxdump.tar.gz` to use a local taxonomy
snapshot instead. If the protein FASTA or taxon map is omitted, PETRA downloads
NCBI NR or the full NCBI protein taxon map. Downloaded archives and DIAMOND
inputs are removed after the build; the extracted taxonomy is retained. The
command fails rather than replacing an existing download or `petra.dmnd`.

Add `--build-rna-dbs` to also build the databases used for rRNA detection, next
to `petra.dmnd` in the same output directory:

- `petra.rrna.cm`, the `cmpress`ed concatenation of the 5S, 5.8S, SSU and LSU
  Rfam families, for `cmscan`.
- `petra.rrna.clanin`, built alongside `petra.rrna.cm` and filtered to the
  families present in that CM, for `cmscan` clan competition.
- `petra.ssu.fasta`, a BLAST database of SILVA SSURef NR99 and PR2, carrying
  NCBI taxon identifiers.
- `petra.lsu.fasta`, a BLAST database of SILVA LSURef NR99 and EUKARYOME,
  carrying NCBI taxon identifiers.

This part requires `makeblastdb`, `cmpress` and `taxonkit` on `PATH`. PETRA
passes the bundled `database/taxonomy/` directory to TaxonKit while constructing
the RNA databases. SILVA sequences of unidentified organisms are discarded, and
PR2 sequences already present in SILVA are dropped so the SSU database holds
each accession once. PR2 and EUKARYOME sequences take the NCBI taxon of the
deepest name in their lineage that TaxonKit recognizes. EUKARYOME
`unclassified` fillers are skipped, and its species epithets are looked up
as binomials under their genus. Other sources are downloaded into the output
directory and removed after the build.

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

To avoid passing the database on each run:

```console
export PETRA_DB_DIR=/path/to/database
petra run --assembly assembly.fasta --long-reads reads.fastq.gz --output run-1
```

#### Taxonomic assignment

PETRA first checks the number of aligned unigenes against a
length-dependent minimum, once for the whole sequence. With the defaults,
sequences up to 50,000 bp require three aligned unigenes. Above that length,
the minimum is `3 + ceil(3 * log10(sequence_length / 50000))`:

| Sequence length | Minimum aligned unigenes |
| --- | --- |
| Up to 50 kb | 3 |
| 100 kb | 4 |
| 500 kb | 6 |
| 5 Mb | 9 |
| 50 Mb | 12 |

Configure the curve with:
- `--min-aligned-unigenes` (default `3`, minimum `3`): the minimum at or below
  the length anchor.
- `--unigene-support-length` (default `50000`, minimum `1`): the anchor in bp.
- `--unigenes-per-decade` (default `3`, minimum `1`): additional unigenes per
  tenfold increase above the anchor.

Sequences with taxonomic evidence below this minimum have all ranks set to
`N/A` and `unassigned_reason` set to `insufficient_support`. Their gene-level
evidence remains in `predictions.gene_taxonomy.tsv`.

Sequences passing the gate assign a taxon only when it holds more than the
`--alpha` share (default `0.5`, a strict majority) of all aligned unigenes.
Ranks are resolved top-down: at each rank, only unigenes that agree with the
consensus at every higher rank vote. Ranks that the winning clade lacks (for
example, the kingdom of Ciliophora) do not stop the descent. If no rank
exceeds `--alpha` despite taxonomic evidence, the reason is `no_consensus`.

#### Optional rRNA and tRNA analyses

Add `--rrna` to detect rRNA genes with `cmscan` and classify SSU and LSU genes
with `blastn`. It requires Infernal's `cmscan` and BLAST's `blastn` on `PATH`
and produces these files under the output directory's `rna/` subdirectory:

- `rna/predictions.rrna.gff`, the full rRNA gene inventory, including 5S and
  5.8S.
- `rna/predictions.rrna.tsv`, the best BLAST hit per SSU or LSU gene, with a
  pipe-joined lineage and a `non_host` column.
- `rna/predictions.rrna.warnings.txt`, host-classification warnings, empty when
  there are none.

`--host-taxon` accepts an NCBI scientific name or TaxID. PETRA resolves the
taxon to its NCBI phylum and matches that phylum as a whole lineage field.
Without a host, `non_host` is `-`. If the host has no phylum, the analysis
continues with `non_host=-` for all hits and records a warning identifying the
host in `rna/predictions.rrna.warnings.txt`. This file is regenerated on reruns.
Unknown or ambiguous host identifiers and TaxonKit failures still stop the run.
When a host phylum is available, a hit whose lineage does not resolve to that
phylum is reported as non-host. Use `--rrna-min-alignment-length` to change
the minimum reported BLAST alignment length from its default of 200.

Add `--trna` to annotate tRNAs. It requires `tRNAscan-SE` on `PATH` and writes
`rna/predictions.trna.tsv` without pseudogenes, undetermined tRNAs or suppressor
tRNAs, plus `rna/predictions.trna.stats`.

Use `--skip-taxonomic-assignment` to skip the protein-based taxonomic
assignment. In this mode PETRA runs Whokaryote, writing
`whokaryote.noneuk.txt`, plus any analyses selected with `--rrna` or `--trna`.
The database is still required for rRNA classification, but runs
without `--rrna` do not require one:

```console
petra run \
  --assembly assembly.fasta \
  --database database \
  --rrna \
  --trna \
  --skip-taxonomic-assignment \
  --output run-1
```

# Output files

The assembly sequences split by final PETRA domain are written under the
output directory's `decontaminated/` subdirectory:

- `decontaminated/predictions.eukaryota.fasta`
- `decontaminated/predictions.bacteria.fasta`
- `decontaminated/predictions.viruses.fasta`
- `decontaminated/predictions.archaea.fasta`
- `decontaminated/predictions.unassigned.fasta`

A full run also writes `predictions.organelles.tsv`,
listing each sequence that Tiara classified as an organelle and its second-stage
type (`mitochondrion`, `plastid`, or `unknown`). Only sequences of at least
5,000 bp are considered for either output.

The protein-taxonomy workflow also writes `predictions.gene_taxonomy.tsv`.
Each assembly sequence is an unindented header followed by one tab-indented row
per aligned unigene, containing the representative gene
identifier and the lineage supported by its highest-scoring retained hit or
hits. The lineage is semicolon-separated, and missing ranks are omitted;
a completely unresolved best hit is reported as `N/A`. Sequences
without an aligned gene are retained as headers with no rows:

```text
contig_1
	contig_1_1	Eukaryota;Metazoa;Chordata
	contig_1_2	Bacteria;Pseudomonadota
contig_2
	contig_2_1	Viruses
```

The run writes `predictions.summary.png`, a summary graph to eyeball the result:
GC against depth of coverage with point area scaled by sequence length, both
cumulative-length marginals, the five longest sequences, and the share of
sequences and bases that received an assignment.

The taxonomy legend lists each taxon with its cumulative length and sequence
count (`Rhodophyta  53.5 Mbp · 55`), from the largest to the smallest, with
`other` and unassigned sequences last. When the list does not fit in one
column it is split in two, read row by row, so the first row holds the two
largest taxa.

Customize it with `--summary-font-size`, `--summary-disable-panel` (repeat to
drop `header`, `gc-histogram`, `coverage-histogram`, `taxonomy-legend`,
`longest-sequences`, `length-legend` or `assignment-status`), and either
`--summary-theme` or `--summary-colors`.

A theme is a color file packaged with PETRA in `src/petra/themes/`:

- `superkingdom` (default): one color per superkingdom.
- `atlasea`: one color per eukaryotic phylum commonly met in ATLASea samples,
  plus the Dinophyceae and Phaeophyceae classes, which have no phylum in the NCBI
  taxonomy, so eukaryotic contaminants stand out from the host. Every
  superkingdom keeps its default color: other eukaryotes, including those
  without an assigned phylum, stay red.

`--summary-colors` replaces the theme with your own two-column TSV of taxonomic
rank and color, in the same format as the theme files:

```text
Actinopteri	#1F78B4
Bacteria	green
```

A sequence then takes the color of the deepest assigned rank listed in the
file, assigned sequences matching no rank are drawn as `other`, and unassigned
sequences keep the unassigned color.

### Redrawing the summary graph

`petra-plot` renders the same graph from an existing table without rerunning
the workflow, for instance after filtering it with `petra-filter`:

```console
petra-plot -i run-1/predictions.final.tsv -o summary.png --theme atlasea
```

It takes a theme with `--theme` or your own color file with `--colors`, and
`--disable-panel`, `--font-size`, `--title`, `--width`, `--height` and `--dpi`.

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
  --theme atlasea \
  --font-size 12 \
  --output comparison.png
```

Every panel shares one color file: a theme with `--theme` (`superkingdom` by
default) or your own TSV with `--colors`.

To keep the legend short, a taxon keeps its own color only if it holds more
than 5% of the bases of at least one assembly. Smaller taxa take their
superkingdom color instead, which is the one the theme uses for that
superkingdom when it lists it. Change the threshold with `--min-share`, or
set it to `0` to color every taxon.

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
