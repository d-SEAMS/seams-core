#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Build the videos page from per-slug Odysee claim files.

A slug with one claim line renders a title, one sentence, and an
/$/embed/ player. A missing file, or a file with only comments, renders
nothing. Claim ids are never invented here. The line format, the
resolution, and the caption file are documented in docs/videos/README.md.

Odysee documents autoplay=true and t= on the embed URL. It does not
document a query parameter that hides related videos or autoplay-next,
so the iframe src is the bare /$/embed/ URL.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

CLAIM = re.compile(r"^([^:\s]+):([A-Za-z0-9]+)$")
CHANNEL = "https://odysee.com/@d-SEAMS:0b23224c14d44b2a90deee60b7cc8c54e2331e9d"
# Title and one sentence for each cut that has an upload. The gallery
# prints both. A topic page prints the player only.
CARDS = {
    "c": (
        "C",
        "The C program reads ice Ih and prints the neighbour count, q6, published CHILL+, and the six-ring count.",
    ),
    "cpp": (
        "C++",
        "The C++ program reads ice Ih and prints the neighbour pairs, q6, published CHILL+, and the six-ring count.",
    ),
    "fortran": (
        "Fortran",
        "The Fortran recording runs the ice Ih census, then the Fortran example on the sI cell.",
    ),
    "python": (
        "Python",
        "The Python program reads ice Ih and prints the neighbour count, q6, published CHILL+, and the six-ring count.",
    ),
    "rust": (
        "Rust",
        "The Rust example reads ice Ih and prints the neighbour count, q6, published CHILL+, and the six-ring count.",
    ),
    "julia": (
        "Julia",
        "The Julia recording runs the ice Ih census, then the Julia example on the sI cell.",
    ),
    "r": (
        "R",
        "The R recording runs the ice Ih census, then the R example on the sI cell.",
    ),
    "lua": (
        "Lua and Fennel",
        "Lua and Fennel print the ice Ih census, the cage counts, and the basal-interface CHILL+ labels.",
    ),
    "cli": (
        "The seams command",
        "The seams command reads ice Ih and prints the neighbour count, q6, published CHILL+, and the six-ring count.",
    ),
    "lammps": (
        "LAMMPS",
        "LAMMPS run_hook scores five calls on the cubic cell, then one step on the ice Ih dump.",
    ),
    "gromacs": (
        "GROMACS",
        "gmx-seams reads the trajectory written by gmx mdrun and prints the ice census.",
    ),
    "openmm": (
        "OpenMM",
        "OpenMM scores the ice Ih dump and prints the hexagonal count, the six-ring count, and q6.",
    ),
    "plumed": (
        "PLUMED",
        "The PLUMED driver scores the ice Ih dump and prints the hexagonal count, the six-ring count, and q6.",
    ),
    "insitu": (
        "In situ",
        "The OpenMM how-to and the in-situ call score the ice Ih dump.",
    ),
    "interface": (
        "Basal interface",
        "Published CHILL+ on the basal interface frame prints hexagonal 352 and water 434.",
    ),
    "rings": (
        "Primitive rings",
        "The recording counts primitive rings on exampleTraj, then the six-rings of ice Ih.",
    ),
    "cages": (
        "Cages",
        "The recording counts hexagonal cages on ice Ih, double-diamond cages on the cubic cell, and the sI, sII, and sH cages.",
    ),
    "brine": (
        "Brine front",
        "The brine tutorial prints the ice front on three frames, then the paper frame of 864 waters and 9 ions.",
    ),
    "nucleation": (
        "Nucleation",
        "The recording prints q6bar for frame 1 and frame 148 of the nucleation trajectory.",
    ),
    "migrate": (
        "Migrate from 1.x",
        "The recording runs the 2.x hydrogen-bond, coordination, fingerprint, pair, and domain commands on exampleTraj.",
    ),
}
# Facts that belong in docs/videos/README.md, and strings that must stay
# off the reader pages. The README phrases are the old page sentences.
PAGE_BANNED = (
    "embed.txt",
    "claim-name",
    "1920",
    "captions.vtt",
    "aac_coder",
    "VTT",
    "render_molecules",
    "no player",
    "coming soon",
    "huggingface.co",
    "terminal session",
    "burned",
)
README_REQUIRED = (
    "docs/videos/<slug>/embed.txt",
    "claim-name:full-claim-id",
    "1920 by 1080",
    "captions.vtt",
    "renders no player",
    "-c:a aac -aac_coder fast",
    "Captions are burned into the upload",
)
BINDINGS_BANNED = (
    "embed.txt",
    "claim-name",
    "no player",
    "captions.vtt",
    "coming soon",
)
HTML_BANNED = (
    "huggingface.co",
    "video coming soon",
    "docs/videos/",
    "captions.vtt",
    "Captions are burned",
    "embed.txt",
    "aac_coder",
    "claim-name",
)

# The twenty cuts. rings-cages keeps its shot list on disk and is not a
# gallery video.
CUTS = (
    "brine",
    "c",
    "cages",
    "cli",
    "cpp",
    "fortran",
    "gromacs",
    "insitu",
    "interface",
    "julia",
    "lammps",
    "lua",
    "migrate",
    "nucleation",
    "openmm",
    "plumed",
    "python",
    "r",
    "rings",
    "rust",
)
# Gallery order: languages, the command, the engines, then the tutorials.
READER = (
    "c",
    "cpp",
    "fortran",
    "python",
    "rust",
    "julia",
    "r",
    "lua",
    "cli",
    "lammps",
    "gromacs",
    "openmm",
    "plumed",
    "insitu",
    "interface",
    "rings",
    "cages",
    "brine",
    "nucleation",
    "migrate",
)
NOTES = ("rings-cages",)
# claim-name:claim-id, checked against the d-SEAMS channel.
CLAIMS = {
    "interface": "dseams-interface:74af20635ceeef7945f8e22cb5b02522a865baa6",
    "rings": "dseams-rings:cfb3b6687dd81847ecb005782df8cf5143617be5",
    "cli": "dseams-cli:033f8d1535cfb7bf1e75b4c651e9789d7fa0ab16",
    "c": "dseams-c-api:2f7235016514160a792d2e416ba29c711d989b8b",
    "cpp": "dseams-cpp:cd1c64a69f60f6bfb6e06a14eee583dd761f4a7c",
    "python": "dseams-python:be29b99d4059bf4afe45ad0a9ca71e20bf11beaa",
    "lammps": "dseams-lammps:6902ffe5fcbe4ec864d4ffffbc502b6dda4fafa7",
    "openmm": "dseams-openmm:da6666299d8adc7b00cba6ad57499b5590503e42",
    "plumed": "dseams-plumed:ddd2f9e27e0a126ef3feba3bb6da8378464e57df",
    "insitu": "dseams-insitu:f57605cad84bbfa722e91546584c9df76751b8b5",
    "rust": "dseams-rust:687055b1ac9dfe7596006cb4f2fa466de9492ca5",
    "julia": "dseams-julia:ec237bd28c57b7bdfa2996ccc420c16d9799176f",
    "r": "dseams-r:e0c4871fe0399cf259531b5961c436abd79acadb",
    "nucleation": "dseams-nucleation:965e5772f700e5b1d91f1b069557eb9e87f1c60b",
    "gromacs": "dseams-gromacs:22838fea81b056f214f908414fde65dd6a42f655",
    "fortran": "dseams-fortran:e15a82c2f40ff1fdf10c89ff6c9f78b6b3fd583b",
    "lua": "dseams-lua:9a7ffc26dff9590174a707db848f5f5d6d7a737f",
    "brine": "dseams-brine:124a3a7fa52cdf7aff34ac00cfba6068e240a947",
    "cages": "dseams-cages:1c258cea9384acd2e21eeec509c5eedc050d08db",
    "migrate": "dseams-migrate:ca1381615800e95d41cc8a379d62abd324800e59",
}

# Topic page, relative to the repo, and the recording that belongs on it.
# lua-fennel is the lua slug. chill-interface is the interface slug.
# rings-cages is not a cut.
TOPICS = {
    "cli": (
        "docs/orgmode/reference/cli.org",
        "docs/orgmode/user/howto/classify-cli.org",
    ),
    "c": ("docs/orgmode/user/howto/bindings.org",),
    "cpp": ("docs/orgmode/user/howto/bindings.org",),
    "r": ("docs/orgmode/user/howto/bindings.org",),
    "python": (
        "docs/orgmode/user/howto/bindings.org",
        "docs/orgmode/reference/python-frame.org",
    ),
    "fortran": ("docs/orgmode/user/howto/bindings.org",),
    "rust": ("docs/orgmode/user/howto/bindings.org",),
    "julia": ("docs/orgmode/user/howto/bindings.org",),
    "lua": ("docs/orgmode/user/howto/bindings.org",),
    "lammps": ("docs/orgmode/user/howto/engine-lammps.org",),
    "gromacs": ("docs/orgmode/user/howto/engine-gromacs.org",),
    "openmm": ("docs/orgmode/user/howto/engine-openmm.org",),
    "plumed": ("docs/orgmode/user/howto/engine-plumed.org",),
    "insitu": (
        "docs/orgmode/user/howto/engine-insitu.org",
        "docs/orgmode/user/tutorials/insitu-run.org",
    ),
    "interface": (
        "docs/orgmode/user/tutorials/interface.org",
        "docs/orgmode/user/howto/chill-plus.org",
    ),
    "rings": ("docs/orgmode/user/tutorials/rings.org",),
    "cages": ("docs/orgmode/user/howto/cages-trace.org",),
    "brine": ("docs/orgmode/user/tutorials/brine.org",),
    "nucleation": ("docs/orgmode/user/tutorials/nucleation.org",),
    "migrate": ("docs/orgmode/user/howto/migrate-1x.org",),
}


def fail(message):
    raise SystemExit(message)


def slugs(videos):
    found = []
    for name in CUTS:
        path = videos / name
        if path.is_dir() and (path / "script.md").is_file():
            found.append(name)
    return found


def read_claim(folder):
    path = folder / "embed.txt"
    if not path.is_file():
        return None
    claim = None
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        match = CLAIM.match(line)
        if match is None or len(match.group(2)) < 16:
            fail("embed.txt in {} needs claim-name:full-claim-id".format(folder.name))
        if claim is not None:
            fail("embed.txt in {} has more than one claim".format(folder.name))
        claim = (match.group(1), match.group(2))
    return claim


def xml(text):
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def iframe(title, name, claim_id):
    src = "https://odysee.com/$/embed/{}:{}".format(name, claim_id)
    return (
        '<iframe title="{}" src="{}" width="960" height="540" '
        'allowfullscreen></iframe>'
    ).format(xml(title), src)


def card(slug):
    try:
        return CARDS[slug]
    except KeyError:
        fail("no reader card for {}".format(slug))


def gallery_block(slug, claim):
    """Title, one sentence, and the player. No claim renders nothing."""
    if claim is None:
        return ""
    title, sentence = card(slug)
    lines = [
        '<div id="video-{}">'.format(slug),
        '<h2 id="{}">{}</h2>'.format(slug, xml(title)),
        "<p>{}</p>".format(xml(sentence)),
        iframe(title, claim[0], claim[1]),
        "</div>",
    ]
    return "\n".join("   " + line for line in lines)


def topic_block(slug, claim):
    """The player alone. The topic page already names the subject."""
    if claim is None:
        return ""
    title, _sentence = card(slug)
    return "   " + iframe(title, claim[0], claim[1])


def embed_rst(slug, claim):
    """One player, or an empty fragment when embed.txt has no claim."""
    body = topic_block(slug, claim)
    if not body:
        return ""
    return ".. raw:: html\n\n{}\n".format(body)


def include_line(slug, page=None):
    """The include a topic page uses for this recording.

    Pages under docs/orgmode/user sit one directory deeper than the
    reference pages, so the relative path is not the same on every page.
    """
    if page is None:
        return ".. include:: ../videos/embeds/{}.rst".format(slug)
    dest = Path("docs/orgmode/videos/embeds") / (slug + ".rst")
    rel = os.path.relpath(dest, Path(page).parent)
    return ".. include:: {}".format(rel.replace(os.sep, "/"))


def write_topic_embeds(videos, dest):
    dest.mkdir(parents=True, exist_ok=True)
    for slug in TOPICS:
        claim = read_claim(videos / slug)
        (dest / (slug + ".rst")).write_text(embed_rst(slug, claim))
    return dest


def shown(videos):
    ready = slugs(videos)
    if set(READER) != set(CUTS):
        fail("reader order does not match the cuts")
    missing = [name for name in ready if name not in READER]
    if missing:
        fail("reader order misses {}".format(", ".join(missing)))
    return [name for name in READER if name in set(ready)]


def players_rst(videos):
    bodies = []
    for slug in shown(videos):
        body = gallery_block(slug, read_claim(videos / slug))
        if body:
            bodies.append(body)
    if not bodies:
        return ""
    return ".. raw:: html\n\n" + "\n".join(bodies) + "\n"


def check_reader_copy(page, readme, bindings):
    """The reader pages stay free of build notes. The README keeps them."""
    if CHANNEL not in page:
        fail("videos page is missing the Odysee channel")
    for banned in PAGE_BANNED:
        if banned in page:
            fail("videos page still has maintainer text ({})".format(banned))
    for required in README_REQUIRED:
        if required not in readme:
            fail("videos README is missing ({})".format(required))
    for banned in BINDINGS_BANNED:
        if banned in bindings:
            fail("bindings page still has maintainer text ({})".format(banned))


def check_player_markup(html):
    """A player is an /$/embed/ iframe. Build notes and placeholders are not."""
    for banned in HTML_BANNED:
        if banned in html:
            fail("player html still has maintainer text ({})".format(banned))
    if "odysee.com" in html.replace("odysee.com/$/embed/", ""):
        fail("player html has an Odysee URL that is not /$/embed/")


def write_players(videos, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(players_rst(videos))
    return dest
