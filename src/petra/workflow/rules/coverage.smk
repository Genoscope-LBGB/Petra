_MINIMAP2_PRESETS = {
    "pacbio": "map-hifi",
    "ont": "map-ont",
}
_SORTED_BAM = "aln.sorted.bam"


rule minimap2:
    """Map long reads to the assembly and produce a sorted primary-alignment BAM."""
    input:
        assembly=config["assembly"],
        reads=config["long_reads"],
    output:
        alignment=(_SORTED_BAM if config["no_bam_cleaning"] else temp(_SORTED_BAM)),
    log:
        "logs/minimap2.log",
    benchmark:
        "benchmarks/minimap2.tsv"
    params:
        preset=_MINIMAP2_PRESETS[config["sequencing_technology"]],
    shell:
        "minimap2 -I 24G -t {threads} -ax {params.preset} "
        "{input.assembly:q} {input.reads:q} 2> {log:q} | "
        "samtools view -F 2308 -u - 2>> {log:q} | "
        "samtools sort -m 1g -@ {threads} -o {output.alignment:q} - "
        "2>> {log:q}"


rule bedtools_genomecov:
    """Calculate depth runs across every assembly sequence."""
    input:
        alignment=_SORTED_BAM,
    output:
        coverage=temp("coverage.bedgraph"),
    log:
        "logs/bedtools_genomecov.log",
    benchmark:
        "benchmarks/bedtools_genomecov.tsv"
    shell:
        "bedtools genomecov -ibam {input.alignment:q} -bga "
        "> {output.coverage:q} 2> {log:q}"


rule median_coverage:
    """Calculate the median per-base coverage of every assembly sequence."""
    input:
        coverage="coverage.bedgraph",
    output:
        coverage=temp("predictions.coverage.tsv"),
    log:
        "logs/median_coverage.log",
    benchmark:
        "benchmarks/median_coverage.tsv"
    script:
        "../scripts/calculate_median_coverage.py"
