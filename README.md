# PETRA

PETRA (Predicted protein-enabled Taxonomic Rank Assignment) finds contaminants
in genome assemblies. It predicts genes on each contig, aligns them to a protein
database with DIAMOND, and assigns the contig to a taxon when most of its genes
agree. Contigs are then split by domain so the host genome can be kept apart.

## Install

```console
pip install .
```

`diamond` must be on `PATH`. The rRNA analyses also need Infernal (`cmscan`,
`cmpress`), BLAST (`blastn`, `makeblastdb`) and `taxonkit`, the tRNA analysis
needs `tRNAscan-SE`.

## Quick start

```console
# 1. Build the database once (downloads NCBI NR and taxonomy)
petra build-db --output database --build-rna-dbs

# 2. Run on an assembly
petra run \
  --assembly assembly.fasta \
  --database database \
  --long-reads reads.fastq.gz \
  --output run-1

# 3. Look at the split by domain and taxonomic assignments
ls run-1/decontaminated/
ls run-1/predictions.final.tsv
```

Set `PETRA_DB_DIR` to skip giving `--database` on every run.

## How contigs are assigned

A contig needs enough aligned genes before PETRA will judge it: 3 up to 50 kb,
then a few more for each tenfold increase in length (4 at 100 kb, 9 at 5 Mb).
Contigs below that are left unassigned (`insufficient_support`).

PETRA then votes rank by rank, from domain down. A taxon wins a rank when it
holds more than half of the genes (`--alpha`, default 0.5), and only genes that
agreed on the ranks above take part in the next vote. If nothing wins, the
contig is unassigned (`no_consensus`).

The thresholds are adjustable with `--min-aligned-unigenes`,
`--unigene-support-length` and `--unigenes-per-decade`.

## Outputs

| File | Content |
| --- | --- |
| `decontaminated/predictions.final.fasta` | Eukaryotic and unassigned contigs: what you usually keep |
| `decontaminated/predictions.{eukaryota,bacteria,archaea,viruses,unassigned}.fasta` | Contigs split by domain |
| `predictions.final.tsv` | One row per contig: length, GC, median coverage, lineage and why it is unassigned, if it is |
| `predictions.gene_taxonomy.tsv` | The lineage of each aligned gene, grouped by contig |
| `predictions.organelles.tsv` | Contigs Tiara flags as mitochondrion or plastid |
| `predictions.summary.png` | GC vs coverage plot, coloured by taxon |
| `rna/` | rRNA and tRNA results, with `--rrna` / `--trna` |

Only contigs of at least 5 kb are considered for the FASTA and organelle outputs.

## Optional analyses

- `--rrna` finds rRNA genes and classifies SSU/LSU hits against SILVA, PR2 and
  EUKARYOME. With `--host-taxon "Chelon auratus"`, hits outside the host's
  phylum are flagged as `non_host`.
- `--trna` annotates tRNAs.
- `--skip-taxonomic-assignment` runs only the RNA analyses and Whokaryote.
- `--keep-whokaryote-euk-headers` keeps the names of sequences Whokaryote
  calls eukaryotic in `whokaryote.euk.txt`.
- `--keep-pyrodigal-predictions` keeps Pyrodigal's predicted proteins in
  `predictions.faa`.
- `--technology ont` switches read mapping from PacBio HiFi to ONT.

## Looking at the results

```console
# Redraw the plot, colouring common eukaryotic phyla
petra-plot -i run-1/predictions.final.tsv -o summary.png --theme atlasea

# Filter the table
petra-filter -i run-1/predictions.final.tsv -o kept.tsv \
  "median_coverage<40 and gc_percent>59 and not phylum==Actinopteri"

# Compare several runs side by side
petra-compare \
  --directory run-1 --name "Chelon auratus" \
  --directory run-2 --name "Jania squamata" \
  --output comparison.png
```

Two themes are bundled: `superkingdom` (default) and `atlasea`, which gives
common eukaryotic phyla their own colours so eukaryotic contaminants stand out
from the host. You can also pass your own colours as a two-column TSV
(`Actinopteri	#1F78B4`) with `--colors`.

## Database details

`petra build-db` creates `petra.dmnd` from NCBI NR (or `--fasta` plus
`--taxon-map`) and keeps the NCBI taxonomy under `database/taxonomy/`. The
taxonomy is also embedded in `petra.dmnd`, so DIAMOND's `--taxonlist` and
`--taxon-exclude` work on it.

- `--with-marferret` adds marine eukaryote proteins from
  [MarFERReT](https://github.com/armbrustlab/marferret) v1.1.1
  ([Groussman et al. 2023](https://doi.org/10.1038/s41597-023-02842-4)).
- `--build-rna-dbs` adds the rRNA databases (Rfam models, SILVA + PR2 for SSU,
  SILVA + EUKARYOME for LSU).
- `--skip-diamond-db` rebuilds only the RNA databases and leaves `petra.dmnd`
  as is.
- `--taxdump` uses a local taxonomy snapshot instead of downloading one.
