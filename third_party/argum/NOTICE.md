# argum

`argum.h` is the single-header amalgamation of
[gershnik/argum](https://github.com/gershnik/argum), Eugene Gershnik's
command-line argument parser. The terms are the BSD 3-Clause License,
Copyright (c) 2022, Eugene Gershnik, in `LICENSE` beside `argum.h`.
That is the file the header comment names; it is upstream argum v2.9's
`LICENSE` unchanged. `src/seams_cli.cpp` uses it for `--help`, `--version` and
every other CLI flag through `Argum::Parser` and `Argum::Colorizer`.

## Why this is not a Meson wrap

argum has no entry in the Meson WrapDB and ships no `meson.build`
(upstream builds with CMake only), so it fails the "exists in Meson
WrapDB, or has a meson-buildable release" test this project otherwise
applies to a vendored third-party header (see
`subprojects/cxxopts.wrap`, `subprojects/backward-cpp.wrap` and
`subprojects/rang.wrap` for headers that pass that test).

More importantly, this copy is not pristine upstream argum: four
follow-up commits patched it in place, and wrapping upstream verbatim
would silently drop all four.

* `49eddbd2` ("fix: make signed index conversions explicit") adds
  `static_cast`s so argum compiles clean under this project's
  `-Wconversion -Wsign-conversion -Wshadow -Werror`
  (`.github/workflows/ci_clang_warn.yml`). Upstream argum does not
  build under those flags as of v2.11.
* `b9926795` ("fix(core): keep the x86 kernels on clang-cl and give
  MSVC a scalar path") lowercases `<Windows.h>` to `<windows.h>` for
  MinGW cross-compilation, where the case-sensitive header lookup only
  has the lowercase name.
* `f31a388a` ("fix: keep the 2.13 colours spellings as aliases") adds
  `BasicOption::alias()` / `OptionNames::addHidden()` / `::visible()`
  so `--colour-types` keeps parsing as a hidden alias of
  `--color-types` without appearing in `--help`. `src/seams_cli.cpp`
  calls `.alias("--colour-types")`; upstream argum has no `alias()`.
* `5dd4e626` ("fix(build): pass iswspace a wint_t in the vendored
  argument parser") fixes one more `-Wconversion` warning the same
  way.

Given that, moving the file next to the code that uses it is more
honest than wrapping a fork and carrying a patch to re-apply it forever
(iralib_interf.h follows the same pattern, living in
`src/include/internal/` since it too has no upstream wrap or
meson-buildable release). This sits under `third_party/`
instead of `src/include/internal/` so it stays out of
`.clang-tidy`'s `HeaderFilterRegex` and Doxygen's `INPUT` the way
`third_party/csm` already does: those enforce this project's own style
on its own code, not on a third-party header.

## Provenance

This copy started as upstream argum v2.9 (commit
`13615b4748272ebe4955dada30f2d4438810757b`) and has the four commits
above applied on top, squashed into `upstream-v2.9.patch`. To
reproduce it:

```bash
git clone --branch v2.9 https://github.com/gershnik/argum.git
cd argum
git apply /path/to/third_party/argum/upstream-v2.9.patch
diff single-file/argum.h /path/to/third_party/argum/argum.h  # empty
```

`upstream-v2.9.patch` is a one-commit mbox, so `git am` (instead of
`git clone` + `git apply`) also applies it directly, with no
whitespace warnings either way.

The full source history, including the exact diff against upstream
v2.9, is this repository's own history for this path and the paths it
was moved from (`git log --follow -- third_party/argum/argum.h`).
