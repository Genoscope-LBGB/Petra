_MINIMAP2_PRESETS = {
    "pacbio": "map-hifi",
    "ont": "map-ont",
}


rule minimap2:
    """Map long reads to the assembly and produce a sorted primary-alignment BAM."""
    input:
        assembly=config["assembly"],
        reads=config["long_reads"],
    output:
        alignment="aln.sorted.bam",
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
