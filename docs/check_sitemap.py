#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Each sitemap URL must name a file written under docs/build.

sphinx-sitemap defaults to ``{lang}{version}{link}``. Sphinx sets language
to ``en``, which listed ``https://docs.dseams.info/en/...`` and those paths
404. ``conf.py`` sets ``html_baseurl`` to ``https://docs.dseams.info/`` and
``sitemap_url_scheme`` to ``{link}``.

Usage:
  docs/check_sitemap.py --self-test
  docs/check_sitemap.py
"""

import ast
import re
import sys
import tempfile
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONF = ROOT / "docs" / "source" / "conf.py"
WORKFLOW = ROOT / ".github" / "workflows" / "ci_docs.yml"
NOTE = ROOT / "docs" / "orgmode" / "development" / "index.org"
BUILD = ROOT / "docs" / "build"
HOST = "docs.dseams.info"
DEPLOY_IF = (
    "github.repository == 'd-SEAMS/seams-core' "
    "&& github.event_name == 'push' "
    "&& github.ref == 'refs/heads/main'"
)


def _local(tag):
    if tag.startswith("{"):
        return tag.partition("}")[2]
    return tag


def locate(loc):
    """Return the docs/build-relative path for one sitemap URL."""
    parsed = urllib.parse.urlsplit(loc.strip())
    if (
        parsed.scheme != "https"
        or parsed.hostname != HOST
        or parsed.port is not None
        or parsed.username
        or parsed.password
        or parsed.query
        or parsed.fragment
    ):
        return None, "{} is outside https://{}/".format(loc, HOST)
    raw = urllib.parse.unquote(parsed.path)
    if "\\" in raw or raw.startswith("//") or "//" in raw:
        return None, "{} has a malformed path".format(loc)
    if raw.startswith("/"):
        raw = raw[1:]
    if raw == "" or raw.endswith("/"):
        raw += "index.html"
    rel = Path(raw)
    if rel.is_absolute() or ".." in rel.parts or any(part == "" for part in rel.parts):
        return None, "{} escapes the build tree".format(loc)
    return rel, None


def canonical_path(rel):
    """Return the public route for one file in the publish tree."""
    route = "/" + rel.as_posix()
    if route.endswith("/index.html"):
        route = route[:-10] or "/"
    elif route.endswith(".html"):
        route = route[:-5]
    return route.lower()


def _urls(url_elem):
    found = []
    for child in list(url_elem):
        name = _local(child.tag)
        if name == "loc" and child.text and child.text.strip():
            found.append(child.text.strip())
        elif name == "link":
            href = child.attrib.get("href", "").strip()
            if href:
                found.append(href)
    return found


def _page_locs(path, build, seen):
    resolved = path.resolve()
    if resolved in seen:
        return [], ["{} repeats a sitemap already read".format(path.name)]
    seen.add(resolved)
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as exc:
        return [], ["{} is not sitemap XML ({})".format(path.name, exc)]
    kind = _local(root.tag)
    if kind == "urlset":
        locs = []
        problems = []
        urls = [child for child in list(root) if _local(child.tag) == "url"]
        if not urls:
            return [], ["{} has no URLs".format(path.name)]
        for url in urls:
            found = _urls(url)
            if not found:
                problems.append("{} has a url entry with no loc".format(path.name))
                continue
            locs.extend(found)
        return locs, problems
    if kind == "sitemapindex":
        locs = []
        problems = []
        entries = [child for child in list(root) if _local(child.tag) == "sitemap"]
        if not entries:
            return [], ["{} has no child sitemaps".format(path.name)]
        for entry in entries:
            found = _urls(entry)
            if not found:
                problems.append("{} has a sitemap entry with no loc".format(path.name))
                continue
            for loc in found:
                rel, err = locate(loc)
                if err:
                    problems.append(err)
                    continue
                child = build.joinpath(*rel.parts)
                if not child.is_file():
                    problems.append("{} does not map to {}".format(loc, rel.as_posix()))
                    continue
                nested, nested_problems = _page_locs(child, build, seen)
                locs.extend(nested)
                problems.extend(nested_problems)
        return locs, problems
    return [], ["{} root is {} rather than urlset".format(path.name, kind)]


def mapped_files(build, locs):
    """Pair each URL with the file it must hit. Problems are separate."""
    problems = []
    files = []
    root = build.resolve()
    routes = {}
    for path in build.rglob("*.html"):
        route = canonical_path(path.relative_to(build))
        if route in routes:
            problems.append(
                "{} and {} have the same canonical route {}".format(
                    routes[route].relative_to(build), path.relative_to(build), route
                )
            )
        routes[route] = path
    for loc in locs:
        rel, err = locate(loc)
        if err:
            problems.append(err)
            continue
        parsed = urllib.parse.urlsplit(loc.strip())
        route = urllib.parse.unquote(parsed.path).lower()
        candidate = routes.get(route)
        if candidate is None and rel.suffix != ".html":
            direct = build.joinpath(*rel.parts)
            if direct.is_file():
                candidate = direct
        if candidate is None:
            problems.append("{} is not a canonical build route".format(loc))
            continue
        try:
            candidate.resolve().relative_to(root)
        except ValueError:
            problems.append("{} escapes the build tree".format(loc))
            continue
        files.append(candidate)
    return files, problems


def sitemap_problems(build):
    sitemap = build / "sitemap.xml"
    if not sitemap.is_file():
        label = "docs/build" if build == BUILD else build.name
        return ["{}/sitemap.xml is missing".format(label)]
    locs, problems = _page_locs(sitemap, build, set())
    if not locs and not any("no URL" in item or "no url" in item for item in problems):
        problems.append("sitemap has no URLs")
    _, mapped = mapped_files(build, locs)
    problems.extend(mapped)
    return problems


def cname_problems(build):
    """GitHub Pages reads CNAME from the publish root and would fight Netlify."""
    if not build.is_dir():
        return []
    found = []
    for path in build.rglob("*"):
        if path.is_file() and path.name.lower() == "cname":
            found.append(path.relative_to(build).as_posix())
    if not found:
        return []
    return ["CNAME in the publish tree: " + ", ".join(found)]


def _assigned(tree):
    found = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if isinstance(target, ast.Name):
            found[target.id] = node.value
    return found


def _dict_strings(node):
    out = {}
    if not isinstance(node, ast.Dict):
        return out
    for key, value in zip(node.keys, node.values):
        if (
            isinstance(key, ast.Constant)
            and isinstance(key.value, str)
            and isinstance(value, ast.Constant)
            and isinstance(value.value, str)
        ):
            out[key.value] = value.value
    return out


def _strings(tree):
    found = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            found.append(node.value)
    return found


def conf_problems():
    text = CONF.read_text(encoding="utf-8")
    tree = ast.parse(text)
    assigned = _assigned(tree)
    problems = []
    base = assigned.get("html_baseurl")
    scheme = assigned.get("sitemap_url_scheme")
    context = assigned.get("html_context")
    if not isinstance(base, ast.Constant) or base.value != "https://docs.dseams.info/":
        problems.append("html_baseurl is not https://docs.dseams.info/")
    if not isinstance(scheme, ast.Constant) or scheme.value != "{link}":
        problems.append('sitemap_url_scheme is not "{link}"')
    fields = _dict_strings(context) if context is not None else {}
    expected = {
        "source_type": "github",
        "source_user": "d-SEAMS",
        "source_repo": "seams-core",
        "source_version": "main",
    }
    for key, value in expected.items():
        if fields.get(key) != value:
            problems.append("html_context {} is not {}".format(key, value))
    for value in _strings(tree):
        if "HaoZeke" in value or "NETLIFY" in value:
            problems.append("conf.py still names a fork deploy target")
            break
    return problems


def workflow_problems():
    text = WORKFLOW.read_text(encoding="utf-8")
    problems = []
    if DEPLOY_IF not in text:
        problems.append("Pages deploy is not limited to d-SEAMS/seams-core main")
    if re.search(r"(?m)^\s*cname\s*:", text):
        problems.append("ci_docs.yml sets a custom domain file")
    if re.search(r"(?mi)^\s*disable_nojekyll\s*:\s*['\"]?true", text):
        problems.append("ci_docs.yml disables .nojekyll")
    if re.search(r"(?mi)^\s*enable_jekyll\s*:\s*['\"]?true", text):
        problems.append("ci_docs.yml enables Jekyll")
    for secret in ("NETLIFY_AUTH_TOKEN", "NETLIFY_SITE_ID", "netlify-cli"):
        if secret in text:
            problems.append("ci_docs.yml still references " + secret)
    block = (
        "  build-deploy:\n"
        "    runs-on: ubuntu-latest\n"
        "    permissions:\n"
        "      contents: write\n"
    )
    if block not in text:
        problems.append("build-deploy is missing permissions: contents: write")
    if "python3 docs/check_mermaid.py" not in text:
        problems.append("Mermaid check is missing")
    if "lycheeverse/lychee-action@" not in text:
        problems.append("Lychee check is missing")
    if "python3 docs/check_sitemap.py\n" not in text:
        problems.append("sitemap file check is missing from ci_docs.yml")
    return problems


def _write_urlset(path, locs):
    rows = ["<?xml version='1.0' encoding='utf-8'?>",
            '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">']
    for loc in locs:
        rows.append("<url><loc>{}</loc></url>".format(loc))
    rows.append("</urlset>\n")
    path.write_text("\n".join(rows), encoding="utf-8")


def _case_tree(root):
    build = root / "build"
    (build / "contributing").mkdir(parents=True)
    (build / "index.html").write_text("home", encoding="utf-8")
    (build / "contributing" / "index.html").write_text("page", encoding="utf-8")
    return build


def mapper_cases():
    problems = []
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        build = _case_tree(root)
        good = [
            "https://docs.dseams.info/",
            "https://docs.dseams.info/contributing/",
        ]
        _write_urlset(build / "sitemap.xml", good)
        got = sitemap_problems(build)
        if got:
            problems.append("valid sitemap was rejected: " + "; ".join(got))
        if cname_problems(build):
            problems.append("valid tree reported a CNAME file")
        (build / "CNAME").write_text("docs.dseams.info\n", encoding="utf-8")
        if not cname_problems(build):
            problems.append("CNAME in the publish tree was accepted")
        (build / "CNAME").unlink()

        bad_cases = [
            (["https://docs.dseams.info/en/index.html"], "en/index.html"),
            (["https://example.com/index.html"], "outside"),
            (["https://docs.dseams.info/missing.html"], "missing.html"),
            (["https://docs.dseams.info/../index.html"], "escapes"),
            (["https://docs.dseams.info/index.html"], "canonical"),
            (["https://docs.dseams.info/contributing/index.html"], "canonical"),
        ]
        for locs, needle in bad_cases:
            _write_urlset(build / "sitemap.xml", locs)
            got = sitemap_problems(build)
            if not got or not any(needle in item for item in got):
                problems.append("bad URL {} was not reported ({})".format(
                    needle, "; ".join(got)))

        _write_urlset(build / "sitemap.xml", [])
        got = sitemap_problems(build)
        if not any("no URL" in item or "no url" in item for item in got):
            problems.append("empty sitemap was accepted")

        pages = root / "pages.xml"
        _write_urlset(pages, ["https://docs.dseams.info/"])
        # The index entry names a file in the build, so copy it there.
        (build / "pages.xml").write_text(pages.read_text(encoding="utf-8"), encoding="utf-8")
        index = (
            "<?xml version='1.0' encoding='utf-8'?>\n"
            '<sitemapindex xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n'
            "<sitemap><loc>https://docs.dseams.info/pages.xml</loc></sitemap>\n"
            "</sitemapindex>\n"
        )
        (build / "sitemap.xml").write_text(index, encoding="utf-8")
        got = sitemap_problems(build)
        if got:
            problems.append("sitemap index was rejected: " + "; ".join(got))
    return problems


def self_test():
    problems = []
    problems.extend(conf_problems())
    problems.extend(workflow_problems())
    problems.extend(mapper_cases())
    return problems


def main(argv):
    if argv[1:] == ["--self-test"]:
        problems = self_test()
    elif not argv[1:]:
        problems = []
        problems.extend(conf_problems())
        problems.extend(workflow_problems())
        problems.extend(cname_problems(BUILD))
        problems.extend(sitemap_problems(BUILD))
    else:
        print("usage: docs/check_sitemap.py [--self-test]", file=sys.stderr)
        return 2
    if problems:
        print("\n".join(problems), file=sys.stderr)
        return 1
    print("sitemap check ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
