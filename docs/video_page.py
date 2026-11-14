#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Build the videos page from per-slug Odysee claim files.

docs/videos/<slug>/embed.txt holds one line, claim-name:full-claim-id,
when that recording has been uploaded. A missing file, or a file with
only comments, renders no player. Claim ids are never invented here.

Odysee documents autoplay=true and t= on the embed URL. It does not
document a query parameter that hides related videos or autoplay-next,
so the iframe src is the bare /$/embed/ URL.
"""

from __future__ import annotations

import re
from pathlib import Path

CLAIM = re.compile(r"^([^:\s]+):([A-Za-z0-9]+)$")
BLURB = {
    "c": "C on ice Ih.",
    "cpp": "C++ on ice Ih.",
    "python": "Python on ice Ih.",
    "fortran": "Fortran on ice Ih.",
    "rust": "Rust on ice Ih.",
    "julia": "Julia on ice Ih.",
    "r": "R on ice Ih.",
    "cli": "The seams command on ice Ih.",
    "interface": "CHILL+ on the basal interface.",
    "rings-cages": "Rings and sI, sII cages.",
    "brine": "The brine front.",
    "nucleation": "Nucleation q6bar.",
    "insitu": "LAMMPS, GROMACS, OpenMM, and PLUMED.",
    "lammps": "LAMMPS run_hook.",
    "gromacs": "GROMACS gmx-seams.",
    "plumed": "PLUMED driver.",
    "lua": "Lua and Fennel.",
    "migrate": "Migrate from 1.x.",
    "openmm": "OpenMM on ice Ih.",
    "rings": "Primitive rings on exampleTraj, then ice Ih.",
    "cages": "Ice Ih hexagonal cages, cubic double-diamond cages, then the clathrate cells.",
}

# The nineteen cuts. migrate and rings-cages keep their shot lists on disk
# and are not gallery videos.
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
    "nucleation",
    "openmm",
    "plumed",
    "python",
    "r",
    "rings",
    "rust",
)
NOTES = ("migrate", "rings-cages")
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
}

# Topic page, relative to the repo, and the recording that belongs on it.
# lua-fennel is the lua slug. chill-interface is the interface slug.
# cages is a gallery video with no topic page. migrate and rings-cages are not cuts.
TOPICS = {
    "cli": (
        "docs/orgmode/reference/cli.org",
        "docs/orgmode/howto/classify-cli.org",
    ),
    "cpp": ("docs/orgmode/howto/bindings.org",),
    "python": (
        "docs/orgmode/howto/bindings.org",
        "docs/orgmode/reference/python-frame.org",
    ),
    "fortran": ("docs/orgmode/howto/bindings.org",),
    "rust": ("docs/orgmode/howto/bindings.org",),
    "julia": ("docs/orgmode/howto/bindings.org",),
    "lua": ("docs/orgmode/howto/lua.org",),
    "lammps": ("docs/orgmode/howto/engine-lammps.org",),
    "gromacs": ("docs/orgmode/howto/engine-gromacs.org",),
    "openmm": ("docs/orgmode/howto/engine-openmm.org",),
    "plumed": ("docs/orgmode/howto/engine-plumed.org",),
    "insitu": (
        "docs/orgmode/howto/engine-insitu.org",
        "docs/orgmode/tutorials/insitu-run.org",
    ),
    "interface": (
        "docs/orgmode/tutorials/interface.org",
        "docs/orgmode/howto/chill-plus.org",
    ),
    "rings": ("docs/orgmode/tutorials/rings.org",),
    "cages": ("docs/orgmode/howto/cages-trace.org",),
    "brine": ("docs/orgmode/tutorials/brine.org",),
    "nucleation": ("docs/orgmode/tutorials/nucleation.org",),
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


def iframe(name, claim_id):
    src = "https://odysee.com/$/embed/{}:{}".format(name, claim_id)
    return (
        '<iframe title="{}" src="{}" width="960" height="540" '
        'allowfullscreen></iframe>'
    ).format(name, src)


def block(slug, claim):
    if claim is None:
        return ""
    blurb = BLURB.get(slug, slug)
    transcript = (
        '<p>Transcript: <a href="../_static/captions/{}.vtt">'
        "docs/videos/{}/captions.vtt</a></p>"
    ).format(slug, slug)
    body = (
        iframe(claim[0], claim[1])
        + "\n   <p>Captions are burned into the upload.</p>\n   "
        + transcript
    )
    return "   {}\n   <p>{}</p>".format(body, blurb)


def embed_rst(slug, claim):
    """One player, or an empty page when embed.txt has no claim."""
    body = block(slug, claim)
    if not body:
        return ""
    return ".. raw:: html\n\n{}\n".format(body)


def include_line(slug):
    return ".. include:: ../videos/embeds/{}.rst".format(slug)


def write_topic_embeds(videos, dest):
    dest.mkdir(parents=True, exist_ok=True)
    for slug in TOPICS:
        claim = read_claim(videos / slug)
        (dest / (slug + ".rst")).write_text(embed_rst(slug, claim))
    return dest


def players_rst(videos):
    bodies = []
    for slug in slugs(videos):
        body = block(slug, read_claim(videos / slug))
        if body:
            bodies.append(body)
    if not bodies:
        return ""
    return ".. raw:: html\n\n" + "\n".join(bodies) + "\n"


def write_players(videos, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(players_rst(videos))
    return dest
