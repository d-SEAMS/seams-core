# Releasing 2.14.0

Nothing in this file uploads. Run the commands after the `v2.14.0` tag
exists on `d-SEAMS/seams-core`. The tag commit is the source every
registry builds. A tag on a fork does not publish. The GitHub release
workflow and the wheel workflow run only when
`github.repository == 'd-SEAMS/seams-core'`. The release uses
`GITHUB_TOKEN`. PyPI uses trusted publishing. No other token is added
for this repository. Publish in the order below: Julia's General
registry resolves `libseams_jll` before it will accept Dseams, and
crates.io resolves `dseams-sys` before it will accept `dseams`. The
other registries do not depend on each other.

Build with `-j2`. Do not install into `/usr/lib`.

## 0. Tag

```sh
git checkout <release-commit>
git tag -a v2.14.0 -m "d-SEAMS 2.14.0"
git push https://github.com/d-SEAMS/seams-core.git v2.14.0
```

That push is what starts `.github/workflows/release.yml` and the
publish job in `.github/workflows/wheels.yml`. Both jobs are skipped
unless the repository is `d-SEAMS/seams-core`.

`packaging/yggdrasil/L/libseams/build_tarballs.jl` pins
`v2_14_0_tag_sha`. If `git rev-parse v2.14.0` is not that constant,
replace the constant with the tag commit before the Yggdrasil pull
request. The PyPI sdist hash in `recipe/recipe.yaml` is filled in step 7,
after the sdist below is the one that will be uploaded.

## 1. PyPI (`pydseamslib` 2.14.0)

The distribution name is `pydseamslib`. The import name is `dseams`.
See `bindings/python/README.md`.

Wheels are built from the sdist by `.github/workflows/wheels.yml`
(`cibuildwheel` in `pyproject.toml`: manylinux_2_28 and macOS 14). A
`v*` tag on `d-SEAMS/seams-core` publishes through the `pypi` GitHub
environment with trusted publishing. The job does not run on a fork.
To build the same artifacts and stop before publish:

```sh
python3 -m build --sdist --outdir dist
python3 -m cibuildwheel dist/pydseamslib-2.14.0.tar.gz --output-dir dist
python3 -m twine check dist/pydseamslib-2.14.0.tar.gz dist/*.whl
```

`twine check` does not upload. The upload is the trusted-publishing
job, not a token on the command line.

## 2. crates.io (`dseams-sys`, then `dseams`)

`vendor/` is generated and is not in git. `cargo publish` checks the
packaged crate. Publish `dseams-sys` and wait until crates.io will serve
it, then publish `dseams`. A dry run of `dseams` before that upload
cannot see `dseams-sys` on the index.

```sh
python3 bindings/rust/seams-sys/vendor.py
cd bindings/rust
cargo package -p dseams-sys
cargo publish --dry-run -p dseams-sys
cargo package -p dseams
cargo publish --dry-run -p dseams
```

`cargo publish --dry-run` does not upload and does not need a registry
token. docs.rs builds both crates for `x86_64-unknown-linux-gnu`.
`dseams-sys` returns immediately when `DOCS_RS` is set, so the docs
build does not compile Eigen.

## 3. Yggdrasil, then Julia General

Copy `packaging/yggdrasil/L/libseams` into a Yggdrasil checkout and open
that pull request. When the build is registered, `libseams_jll` 2.14.0
is in General.

Then, from a checkout of the `v2.14.0` tag:

```julia
# in a comment on d-SEAMS/seams-core
# @JuliaRegistrator register subdir=bindings/julia/Dseams
```

`.github/workflows/julia-register.yml` posts that comment with
`GITHUB_TOKEN`, and only when
`github.repository == 'd-SEAMS/seams-core'`. TagBot in
`.github/workflows/julia-tagbot.yml` uses the same guard and the same
token.

`bindings/julia/Dseams/Project.toml` already holds the uuid, version
2.14.0, and `[compat]`, with `libseams_jll = "=2.14.0"`.

Until the JLL is registered, test against a local build:

```sh
julia --project=bindings/julia/Dseams bindings/julia/Dseams/dev/develop_jll.jl
# SEAMS_LANG_LIB names libseams_lang from `julia build_tarballs.jl --deploy=local`
# or from a meson install
julia --project=bindings/julia/Dseams -e 'using Pkg; Pkg.test()'
```

## 4. CRAN (`dseams` 2.14.0)

```sh
cd bindings/r/dseams
sh tools/vendor.sh
R CMD build .
R CMD check --as-cran dseams_2.14.0.tar.gz
```

Submit `dseams_2.14.0.tar.gz` at <https://cran.r-project.org/submit.html>.
`cran-comments.md` records the check. The expected CRAN note on a first
upload is `New submission` and nothing else.

## 5. LuaRocks (`dseams` 2.14.0-1)

The rock links an installed `libseams`.

```sh
export SEAMS_DIR="$prefix" SEAMS_LIBDIR="$prefix/lib"
export LD_LIBRARY_PATH="$prefix/lib"
luarocks make bindings/lua/dseams-2.14.0-1.rockspec \
  SEAMS_DIR="$SEAMS_DIR" SEAMS_LIBDIR="$SEAMS_LIBDIR"
luarocks test bindings/lua/dseams-2.14.0-1.rockspec
luarocks pack dseams 2.14.0-1
```

Run the make, test, and pack once for Lua 5.4 and once for LuaJIT
(`--lua-version=5.4` and `--lua-version=5.1` with the LuaJIT interpreter).
Packing does not upload and does not need a LuaRocks token.

## 6. fpm registry

`fpm.toml` at the repository root is the manifest the registry reads.
`packaging/fpm/registry.toml` is the stanza for
<https://github.com/fortran-lang/fpm-registry> `registry.toml`.

The fpm 0.13 registry upload, from `bindings/fortran` (or the root, which
is the same package):

```sh
export FPM_LDFLAGS="-L$prefix/lib -Wl,-rpath,$prefix/lib"
export LD_LIBRARY_PATH="$prefix/lib"
fpm test
fpm build --dump fpm_model.json
fpm publish --show-upload-data
```

`--show-upload-data` does not upload and does not take a token.
Open the registry.toml pull request as well, after `v2.14.0` exists
on `d-SEAMS/seams-core`, because that index fetches the manifest from
the tag.

## 7. conda-forge

After the PyPI sdist from step 1 is the uploaded file:

```sh
sha256sum dist/pydseamslib-2.14.0.tar.gz
```

Write that digest into `recipe/recipe.yaml` (`source.sha256`). The url is
already the PyPI sdist. Then:

```sh
rattler-build build --recipe recipe -c conda-forge
```

Copy `recipe/` into a conda-forge/staged-recipes pull request. The recipe
builds `libseams` and `pydseamslib`. `rattler-build build` is the local
proof. Opening the pull request is the publish step.

## 8. rOpenSci

`docs/ropensci-submission.md` is the software-review issue text. Open it
by hand on <https://github.com/ropensci/software-review/issues/new> after
`pkgcheck::pkgcheck()` is clean. Do not open it as part of tagging.

rOpenSci can co-submit an accepted package to the Journal of Open Source
Software when the package has a short paper. `paper/` is that paper. Say
so in the review thread if you want the co-submission. It is optional and
separate from the CRAN upload.
