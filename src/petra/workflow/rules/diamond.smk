rule diamond:
    """Align clustered predicted proteins against the configured DIAMOND database."""
    input:
        proteins="predictions.clustered.faa",
        database=config["diamond_db"],
    output:
        alignment=temp("predictions.diamond.tsv"),
    log:
        "logs/diamond.log",
    benchmark:
        "benchmarks/diamond.tsv"
    shell:
        "diamond blastp --query {input.proteins:q} --db {input.database:q} "
        "--out {output.alignment:q} --threads {threads} -k 10 "
        "--outfmt 6 qseqid sseqid evalue bitscore staxids "
        "> {log:q} 2>&1"


rule filter_overrepresented_targets:
    """Remove alignments to target sequences occurring above the configured limit."""
    input:
        alignment="predictions.diamond.tsv",
    output:
        alignment=temp("predictions.diamond.filtered.tsv"),
    log:
        "logs/filter_overrepresented_targets.log",
    benchmark:
        "benchmarks/filter_overrepresented_targets.tsv"
    params:
        maximum=config["max_target_occurrences"],
    script:
        "../scripts/filter_overrepresented_targets.py"


rule replace_taxids_with_lineages:
    """Replace subject TaxIDs with one fixed-rank lineage.

    DIAMOND can report multiple, numerically ordered TaxIDs for one subject.
    Their LCA is the most specific lineage supported by every assignment and
    avoids treating the first TaxID as a preferred assignment.
    """
    input:
        alignment="predictions.diamond.filtered.tsv",
    output:
        alignment=temp("predictions.diamond.lineages.tsv"),
    params:
        taxonomy=config["taxonomy"],
    log:
        "logs/replace_taxids_with_lineages.log",
    benchmark:
        "benchmarks/replace_taxids_with_lineages.tsv"
    shell:
        "(taxonkit --data-dir {params.taxonomy:q} --threads {threads} "
        "lca --taxids-field 5 --separator ';' "
        "{input.alignment:q} "
        "| taxonkit --data-dir {params.taxonomy:q} --threads {threads} "
        "reformat2 --taxid-field 6 "
        "--format '{{domain|acellular root|superkingdom}};{{kingdom}};"
        "{{phylum}};{{subphylum}};{{superclass}};{{class}};{{superorder}};"
        "{{order}};{{suborder}};{{family}};{{subfamily}};{{genus}};"
        "{{subgenus}};{{species}};{{subspecies}}' "
        "--miss-rank-repl 'N/A' --miss-taxid-repl 'N/A' "
        "| cut -f 1-4,7) > {output.alignment:q} 2> {log:q}"


rule assign_taxonomy:
    """
    Require length-dependent cluster support once per sequence, then assign its
    taxonomy by the observed alpha-majority of its aligned clusters.
    """
    input:
        assembly=config["assembly"],
        proteins="predictions.clustered.faa",
        alignment="predictions.diamond.lineages.tsv",
    output:
        taxonomy=temp("predictions.taxonomy.raw.tsv"),
        reasons=temp("predictions.unassigned_reasons.assign_taxonomy.tsv"),
        gene_taxonomy="predictions.gene_taxonomy.tsv",
    log:
        "logs/assign_taxonomy.log",
    benchmark:
        "benchmarks/assign_taxonomy.tsv"
    params:
        alpha=config["alpha"],
        min_aligned_clusters=config["min_aligned_clusters"],
        cluster_support_length=config["cluster_support_length"],
        clusters_per_decade=config["clusters_per_decade"],
    script:
        "../scripts/assign_taxonomy.py"
