#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Build the videos page from per-slug Odysee claim files.

docs/videos/<slug>/embed.txt holds one line, claim-name:full-claim-id,
when that recording has been uploaded. A file with only comments means
the upload is not available yet. Claim ids are never invented here.
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
    "rings": "Primitive rings on the mW cell.",
    "cages": "sI, sII, and sH cages.",
}

# Topic page, relative to the repo, and the recording that belongs on it.
# lua-fennel is the lua slug. chill-interface is the interface slug.
# c, r, and migrate have no recording yet.
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
    for path in sorted(videos.iterdir()):
        if path.is_dir() and (path / "script.md").is_file():
            found.append(path.name)
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
    blurb = BLURB.get(slug, slug)
    transcript = (
        '<p>Transcript: <a href="../_static/captions/{}.vtt">'
        "docs/videos/{}/captions.vtt</a></p>"
    ).format(slug, slug)
    if claim is None:
        body = "<p><strong>video coming soon</strong></p>\n   " + transcript
    else:
        body = (
            iframe(claim[0], claim[1])
            + "\n   <p>Captions are burned into the upload.</p>\n   "
            + transcript
        )
    return "   {}\n   <p>{}</p>".format(body, blurb)


def embed_rst(slug, claim):
    """One player, the same claim rule as the gallery."""
    return ".. raw:: html\n\n{}\n".format(block(slug, claim))


def include_line(slug):
    return ".. include:: ../videos/embeds/{}.rst".format(slug)


def write_topic_embeds(videos, dest):
    dest.mkdir(parents=True, exist_ok=True)
    for slug in TOPICS:
        claim = read_claim(videos / slug)
        (dest / (slug + ".rst")).write_text(embed_rst(slug, claim))
    return dest


def players_rst(videos):
    parts = [".. raw:: html", ""]
    for slug in slugs(videos):
        parts.append(block(slug, read_claim(videos / slug)))
    parts.append("")
    return "\n".join(parts)


def write_players(videos, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(players_rst(videos))
    return dest
