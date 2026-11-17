#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""End cards match, and docs iframes stay on the embed player."""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
END = """## End card

Goswami, R., Goswami, A., & Singh, J. K. (2020). d-SEAMS: Deferred Structural Elucidation Analysis for Molecular Simulations. Journal of Chemical Information and Modeling. https://doi.org/10.1021/acs.jcim.0c00031

d-SEAMS 2.0: Incremental ring and cage algorithms for ice structure analysis. Rohit Goswami, Amrita Goswami, Ruhila Goswami, Jayant K. Singh, to be submitted soon.

Questions, issues, discussions and pull requests welcome
"""
IFRAME = re.compile(r"<iframe\b[^>]*\bsrc=[\"']([^\"']+)[\"']", re.I)
LINK = re.compile(r"https?://[^\s)\]>\"']+")
CHANNEL = "https://odysee.com/@d-SEAMS:0b23224c14d44b2a90deee60b7cc8c54e2331e9d"
EMBED = re.compile(r"^https://odysee\.com/\$/embed/[^/\s:]+:[A-Za-z0-9]+$")
SELF = Path(__file__).resolve()


def fail(message):
    print(message, file=sys.stderr)
    raise SystemExit(1)


def pages():
    found = [ROOT / "README.md"]
    docs = ROOT / "docs"
    for path in docs.rglob("*"):
        if not path.is_file() or path.resolve() == SELF:
            continue
        if path.suffix in {".md", ".org", ".rst", ".py"}:
            found.append(path)
    return found


def main():
    scripts = sorted((ROOT / "docs" / "videos").glob("*/script.md"))
    if len(scripts) < 19:
        fail("expected every video script, found {}".format(len(scripts)))
    for path in scripts:
        text = path.read_text(encoding="utf-8")
        if not text.endswith(END):
            fail("{} end card does not match".format(path.relative_to(ROOT)))
        if text.count("## End card") != 1:
            fail("{} end card is repeated".format(path.relative_to(ROOT)))
    for path in pages():
        text = path.read_text(encoding="utf-8", errors="replace")
        for src in IFRAME.findall(text):
            if "odysee.com/@" in src or re.search(r"odysee\.com/(?!\$/embed/)", src):
                fail("{} iframe is a watch URL: {}".format(path.relative_to(ROOT), src))
        for url in LINK.findall(text):
            if "rgvids" in url.lower():
                fail("{} links to @rgvids: {}".format(path.relative_to(ROOT), url))
            if "odysee.com" not in url:
                continue
            bare = url.rstrip(".,")
            embed_prefix = CHANNEL[:CHANNEL.index("@")] + "$/embed/"
            if bare == CHANNEL or EMBED.match(bare):
                continue
            # Builders and the README write the embed form with a placeholder.
            if bare.startswith(embed_prefix) and ("{" in bare or "<" in bare):
                continue
            fail("{} has an Odysee URL that is not the channel or an embed: {}".format(
                path.relative_to(ROOT), url))
    print("video cards ok ({} scripts)".format(len(scripts)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
