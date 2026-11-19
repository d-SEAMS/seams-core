# Contributing to seams-core

This file is the working contract for anyone changing this repository. The
README is for users of the engine; this file is for people editing it.
`docs/orgmode/development/index.org` is the long form, and
`docs/orgmode/development/repositories.org` maps the other d-SEAMS
repositories.

## What lives here and what does not

seams-core holds the C++ engine (`libyodaLib`), the stable C API
`libseams` (`include/seams.h`), the C ABI `libseams_lang` that the
language bindings link (`bindings/cxx/seams_lang.h`), the `seams` CLI, the
LAMMPS, GROMACS, OpenMM, PLUMED and in-situ engine drivers, the C and C++
documentation, the science scripts, the conda recipe and the Yggdrasil
recipe for `libseams_jll`.

- A new analysis, a changed result, a reader for a new file format, a C
  API or `seams_lang.h` entry point, or an engine driver: change it here,
  with a C test and an entry on `docs/api/c.md` or `docs/api/cpp.md`.
- A wrapper for an existing C entry point in Python, Lua or Fennel, R,
  Julia, Rust or Fortran belongs in that binding's repository:
  d-SEAMS/PydSEAMSlib, d-SEAMS/yodaStruct, d-SEAMS/dseams-r,
  d-SEAMS/Dseams.jl, d-SEAMS/dseams-rs, d-SEAMS/dseams-fortran. Land the
  engine change here first; the binding then moves its pin.
- Do not add a language binding, a binding package manifest, or Rust,
  Julia, Lua, Fortran or R code to the docs pages here. A docs page that
  needs one links the binding repository.

## Layout

| Path | Contents |
|---|---|
| `src/` | The engine (`libyodaLib`), the C API implementation and the CLI |
| `src/include/internal/` | Engine headers; not installed |
| `include/` | Installed headers: `seams.h`, `seams.hpp`, `seams_lammps_box.h` |
| `bindings/cxx/` | `seams_lang.h` and `libseams_lang`, the C ABI the bindings link |
| `engines/`, `adapters/python/` | Engine drivers and the ctypes-only `dseams_insitu` reporter |
| `tests/` | Catch2 suites, C API tests, CLI and walk checks, fuzzers |
| `docs/orgmode/` | The book, in Org; exported to `docs/source/` before Sphinx |
| `docs/api/`, `docs/engines/`, `docs/tutorials/` | Markdown pages copied into the book |
| `docs/examples/` | C, C++, CLI and ctypes examples the docs tests run |
| `docs/videos/` | Recording tapes, captions and Odysee embeds |
| `changelog.d/` | Unreleased changelog fragments |
| `recipe/` | conda-forge `libseams` recipe (rattler-build) |
| `packaging/yggdrasil/` | Yggdrasil recipe for `libseams_jll` |
| `subprojects/` | Meson wraps (linkcell, minimage, readcon-core, readcon-db, vesin, ...) |

## Build and test

pixi gives the locked toolchain; Meson builds.

```bash
pixi run setup        # meson setup bbdir -Dwith_tests=true
pixi run build        # ninja -C bbdir
pixi run test         # meson test -C bbdir
```

Or by hand inside `pixi shell` or `nix develop`:

```bash
meson setup bbdir -Dwith_tests=true
meson compile -C bbdir
meson test -C bbdir --print-errorlogs
meson test -C bbdir --suite docs          # documentation checks only
meson test -C bbdir capi_suite capi_abi   # one or more tests by name
```

Tests to know by name:

| Test | Checks |
|---|---|
| `capi_suite`, `capi_abi`, `capi_contract`, `capi_symbols` | The C API, its ABI and its symbol list |
| `seams_lang_fixtures`, `seams_lang_frames` | The `libseams_lang` C ABI the bindings rely on |
| `docs_snippet_<page>` | Every C, C++ and Python block on that docs page |
| `docs_check`, `docs_api_coverage`, `docs_c_symbols` | Every C and C++ symbol is documented with an example |
| `release_versions`, `release_versions_selftest` | Release metadata and its parser |
| `docs_llms_txt` | `llms.txt` matches `scripts/gen_llms_txt.py` |
| `docs_video_tapes` | The recording tapes whose sources live here replay and match `recorded.txt` |

The Python snippets in the docs import `dseams` and `pydseams` from the
installed `pydseamslib` of the same version
(`pip install pydseamslib`, or a build of d-SEAMS/PydSEAMSlib).
A `docs_snippet_` test fails, naming that pip line, when the package is
missing. A test that exits 77 is a skip because an optional toolchain or
input is absent.

Do not weaken, skip or delete a test to make the suite pass. Fix the code
under test. A test is removed only when its subject moved to a binding
repository, and the commit says which.

## For coding agents

`skills/d-seams/SKILL.md` is the user-facing skill: install, the `seams`
command, Python, the C API, and the usual mistakes. `llms.txt` is the
index of the published manual. The docs build also writes `llms-full.txt`
and a `.md` file beside each HTML page. `context7.json` names the folders
Context7 should index.

Read this file before editing the repository. Do not weaken, skip, or
delete a test to make the suite pass.

## Documentation

Edit the Org sources under `docs/orgmode/` and the Markdown pages under
`docs/api/`, `docs/engines/` and `docs/tutorials/`. `docs/source/*.rst`
copies of those pages are generated. `pixi run -e docs docbld` builds the
site. After adding, renaming or retitling a page, regenerate `llms.txt`:

```bash
python3 scripts/gen_llms_txt.py
```

## Commits and changelog

- Conventional Commits: `feat:`, `fix:`, `docs:`, `test:`, `build:`,
  `ci:`, `refactor:`, `perf:`, `chore:`, with a scope where it helps
  (`fix(cages): ...`). One logical change per commit.
- Every user-visible change adds a fragment `changelog.d/+slug.<type>.md`
  with type `added`, `changed`, `deprecated`, `removed`, `fixed` or
  `security`, in the same commit. The release cut renders them with
  towncrier.
- ASCII punctuation in source, docs and commit messages.

## Releases

`docs/orgmode/contributing/developer/release.org`.
