_SUMMARY_COLORS = [config["summary_colors"]] if config.get("summary_colors") else []


rule plot_summary_graph:
    """Render the summary graph used to eyeball a finished run."""
    input:
        taxonomy="predictions.final.tsv",
        colors=_SUMMARY_COLORS,
    output:
        graph="predictions.summary.png",
    log:
        "logs/plot_summary_graph.log",
    benchmark:
        "benchmarks/plot_summary_graph.tsv"
    params:
        font_size=config.get("summary_font_size", 10.0),
        disabled_panels=config.get("summary_disabled_panels", []),
    script:
        "../scripts/plot_summary_graph.py"
