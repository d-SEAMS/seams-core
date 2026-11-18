# Copyright (c) d-SEAMS developers
"""Publish a refresh page for every docname the book has moved.

``page_redirects`` maps the old docname to the new one. The HTML is a
short page with a refresh, so an old URL still opens the page that
replaced it. The stub is not indexed and its canonical link names the
stub itself, so a crawler never treats it as a second copy of the target.
"""

import posixpath


def _target(old, new):
    base = posixpath.dirname(old + ".html") or "."
    return posixpath.relpath(new + ".html", start=base)


def collect_redirects(app):
    mapping = app.config.page_redirects or {}
    base = app.config.html_baseurl or ""
    for old in sorted(mapping):
        new = mapping[old]
        yield (
            old,
            {
                "redirect_target": _target(old, new),
                "redirect_url": base + old + ".html",
                "redirect_title": "Moved",
            },
            "redirect.html",
        )


def setup(app):
    app.add_config_value("page_redirects", {}, "html")
    app.connect("html-collect-pages", collect_redirects)
    return {"version": "1.0", "parallel_read_safe": True, "parallel_write_safe": True}
