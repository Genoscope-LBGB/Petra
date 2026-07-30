from pathlib import Path


rule hello:
    output:
        "results/message.txt"
    run:
        destination = Path(output[0])
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(f"{config['message']}\n", encoding="utf-8")
