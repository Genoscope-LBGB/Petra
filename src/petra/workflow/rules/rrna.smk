_RRNA_CHUNK_INDICES = list(range(config["rrna_chunks"]))
_RRNA_CHUNKS = expand("rrna/chunk.{index}.fasta", index=_RRNA_CHUNK_INDICES)
_RRNA_TBLOUTS = expand("rrna/chunk.{index}.tblout", index=_RRNA_CHUNK_INDICES)
# The rank format is spelled out in the shell command, as in rules/diamond.smk,
# because Snakemake reads '{rank}' in params as a wildcard.


rule split_assembly_for_cmscan:
    """Split the assembly into base-balanced chunks to scatter cmscan."""
    input:
        assembly=config["assembly"],
    output:
        chunks=temp(_RRNA_CHUNKS),
    log:
        "logs/split_assembly_for_cmscan.log",
    benchmark:
        "benchmarks/split_assembly_for_cmscan.tsv"
    script:
        "../scripts/split_assembly_for_cmscan.py"


rule cmscan:
    """Search one chunk for the rRNA families of the covariance model.

    --clanin lets cmscan mark the overlapping hits of the per-domain SSU and LSU
    models as clan competitors, so one rRNA region is not reported several times.
    An empty chunk occurs when the assembly holds fewer sequences than chunks.
    """
    input:
        chunk="rrna/chunk.{index}.fasta",
        model=config["covariance_model"],
        clan=config["clan_info"],
    output:
        tblout=temp("rrna/chunk.{index}.tblout"),
    log:
        "logs/cmscan.{index}.log",
    benchmark:
        "benchmarks/cmscan.{index}.tsv"
    shell:
        "if [ -s {input.chunk:q} ]; then "
        "cmscan --cpu {threads} --cut_ga --rfam --nohmmonly "
        "--clanin {input.clan:q} --fmt 2 --tblout {output.tblout:q} "
        "-o /dev/null {input.model:q} {input.chunk:q} > {log:q} 2>&1; "
        "else : > {output.tblout:q} && : > {log:q}; fi"


rule merge_cmscan_tblout:
    """Concatenate the scattered cmscan tables."""
    input:
        tblout=_RRNA_TBLOUTS,
    output:
        tblout=temp("predictions.rrna.tblout"),
    log:
        "logs/merge_cmscan_tblout.log",
    benchmark:
        "benchmarks/merge_cmscan_tblout.tsv"
    shell:
        "cat {input.tblout:q} > {output.tblout:q} 2> {log:q}"


rule extract_rrna_genes:
    """Write the rRNA gene inventory and the SSU and LSU sequences to classify."""
    input:
        assembly=config["assembly"],
        tblout="predictions.rrna.tblout",
    output:
        gff="predictions.rrna.gff",
        genes=temp("predictions.rrna.genes.tsv"),
        ssu=temp("predictions.rrna.ssu.fasta"),
        lsu=temp("predictions.rrna.lsu.fasta"),
    log:
        "logs/extract_rrna_genes.log",
    benchmark:
        "benchmarks/extract_rrna_genes.tsv"
    script:
        "../scripts/extract_rrna_genes.py"


rule blastn_ssu:
    """Align the small subunit genes against the SILVA and PR2 database."""
    input:
        query="predictions.rrna.ssu.fasta",
        database=config["ssu_database"],
    output:
        alignment=temp("predictions.rrna.ssu.blast.tsv"),
    log:
        "logs/blastn_ssu.log",
    benchmark:
        "benchmarks/blastn_ssu.tsv"
    shell:
        "if [ -s {input.query:q} ]; then "
        "blastn -query {input.query:q} -db {input.database:q} "
        "-out {output.alignment:q} -num_threads {threads} "
        '-outfmt "6 std qlen slen qcovs staxid" > {log:q} 2>&1; '
        "else : > {output.alignment:q} && : > {log:q}; fi"


rule blastn_lsu:
    """Align the large subunit genes against the SILVA and EUKARYOME database."""
    input:
        query="predictions.rrna.lsu.fasta",
        database=config["lsu_database"],
    output:
        alignment=temp("predictions.rrna.lsu.blast.tsv"),
    log:
        "logs/blastn_lsu.log",
    benchmark:
        "benchmarks/blastn_lsu.tsv"
    shell:
        "if [ -s {input.query:q} ]; then "
        "blastn -query {input.query:q} -db {input.database:q} "
        "-out {output.alignment:q} -num_threads {threads} "
        '-outfmt "6 std qlen slen qcovs staxid" > {log:q} 2>&1; '
        "else : > {output.alignment:q} && : > {log:q}; fi"


rule select_rrna_best_hits:
    """Keep one hit per rRNA gene, breaking bitscore ties deterministically."""
    input:
        genes="predictions.rrna.genes.tsv",
        ssu="predictions.rrna.ssu.blast.tsv",
        lsu="predictions.rrna.lsu.blast.tsv",
    output:
        hits=temp("predictions.rrna.hits.tsv"),
    log:
        "logs/select_rrna_best_hits.log",
    benchmark:
        "benchmarks/select_rrna_best_hits.tsv"
    params:
        minimum=config["rrna_min_alignment_length"],
    script:
        "../scripts/select_rrna_best_hits.py"


rule replace_rrna_taxids_with_lineages:
    """Replace the subject TaxID of every retained hit with a fixed-rank lineage."""
    input:
        hits="predictions.rrna.hits.tsv",
    output:
        lineages=temp("predictions.rrna.lineages.tsv"),
    log:
        "logs/replace_rrna_taxids_with_lineages.log",
    benchmark:
        "benchmarks/replace_rrna_taxids_with_lineages.tsv"
    shell:
        "if [ -s {input.hits:q} ]; then "
        "taxonkit --threads {threads} reformat2 --taxid-field 13 "
        "--format '{{domain|acellular root|superkingdom}};{{kingdom}};"
        "{{phylum}};{{subphylum}};{{superclass}};{{class}};{{superorder}};"
        "{{order}};{{suborder}};{{family}};{{subfamily}};{{genus}};"
        "{{subgenus}};{{species}};{{subspecies}}' "
        "--miss-rank-repl 'N/A' --miss-taxid-repl 'N/A' {input.hits:q} "
        "> {output.lineages:q} 2> {log:q}; "
        "else : > {output.lineages:q} && : > {log:q}; fi"


rule summarize_rrna_hits:
    """Report the classified rRNA genes and flag the ones outside the host."""
    input:
        lineages="predictions.rrna.lineages.tsv",
    output:
        taxonomy="predictions.rrna.tsv",
    log:
        "logs/summarize_rrna_hits.log",
    benchmark:
        "benchmarks/summarize_rrna_hits.tsv"
    params:
        host_phylum=config["host_phylum"],
    script:
        "../scripts/summarize_rrna_hits.py"
