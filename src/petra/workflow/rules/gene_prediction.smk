rule pyrodigal:
    """
    -m masks runs of N so genes cannot span them
    -c prevents genes from running off sequence edges
    -p meta uses pre-trained models from different organisms instead of training on the input

    """
    input:
        assembly=config["assembly"],
    output:
        gff=temp("predictions.gff"),
        proteins=temp("predictions.faa"),
    log:
        "logs/pyrodigal.log",
    benchmark:
        "benchmarks/pyrodigal.tsv"
    shell:
        "pyrodigal -m -c -f gff -p meta -j {threads} "
        "-i {input.assembly:q} -o {output.gff:q} -a {output.proteins:q} "
        "> {log:q} 2>&1"


rule segmasker:
    """Mask low-complexity regions in the predicted protein sequences."""
    input:
        proteins="predictions.faa",
    output:
        proteins=temp("predictions.masked.faa"),
    log:
        "logs/segmasker.log",
    benchmark:
        "benchmarks/segmasker.tsv"
    shell:
        "segmasker -in {input.proteins:q} -out {output.proteins:q} "
        "-outfmt fasta > {log:q} 2>&1"


rule filter_low_complexity:
    """Discard proteins exceeding the configured low-complexity threshold."""
    input:
        proteins="predictions.masked.faa",
    output:
        proteins=temp("predictions.filtered.faa"),
    log:
        "logs/filter_low_complexity.log",
    benchmark:
        "benchmarks/filter_low_complexity.tsv"
    params:
        maximum=config["max_low_complexity"],
    script:
        "../scripts/filter_low_complexity.py"


rule cdhit:
    """
    Cluster proteins independently within each originating assembly sequence and
    output at most ``max_clusters`` per assembly sequence.

    """
    input:
        proteins="predictions.filtered.faa",
    output:
        proteins=temp("predictions.clustered.faa"),
        cluster_sizes=temp("predictions.cluster_sizes.tsv"),
    log:
        "logs/cdhit.log",
    benchmark:
        "benchmarks/cdhit.tsv"
    params:
        maximum=config["max_clusters_per_sequence"],
    script:
        "../scripts/cluster_proteins.py"
