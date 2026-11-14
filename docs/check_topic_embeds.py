#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Topic pages use the same embed.txt claim as the gallery.

An iframe src is https://odysee.com/$/embed/<claim-name>:<claim-id>.
A slug with no claim renders nothing. Claim ids are not invented.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from video_page import TOPICS, embed_rst, include_line, read_claim  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
VIDEOS = ROOT / "docs" / "videos"
COMING = ("cages",)


def fail(message):
    print(message, file=sys.stderr)
    raise SystemExit(1)


def main():
    if (VIDEOS / "lua-fennel").exists():
        fail("lua-fennel is the lua slug, not its own directory")
    for slug, pages in TOPICS.items():
        folder = VIDEOS / slug
        if not (folder / "script.md").is_file():
            fail("missing video slug {}".format(slug))
        claim = read_claim(folder)
        html = embed_rst(slug, claim)
        if "odysee.com/@" in html or "/watch" in html:
            fail("{} uses a watch URL".format(slug))
        if claim is None:
            if html or "video coming soon" in html or "odysee.com" in html:
                fail("{} should render nothing".format(slug))
        else:
            src = "https://odysee.com/$/embed/{}:{}".format(claim[0], claim[1])
            if '/$/embed/' not in src or 'src="{}"'.format(src) not in html:
                fail("{} iframe src does not match embed.txt".format(slug))
            if "odysee.com/" in html.replace("odysee.com/$/embed/", ""):
                fail("{} has an Odysee URL that is not /$/embed/".format(slug))
        needle = include_line(slug)
        for page in pages:
            text = (ROOT / page).read_text(encoding="utf-8")
            if needle not in text:
                fail("{} does not embed {}".format(page, slug))
    for slug in COMING:
        if read_claim(VIDEOS / slug) is not None:
            fail("{} has a claim but no recording was supplied".format(slug))
    print("topic embeds ok ({} slugs)".format(len(TOPICS)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
