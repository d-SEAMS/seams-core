# Releasing 2.14.0

This file is the seams-core release: the tag, the GitHub release, conda
`libseams`, and `libseams_jll`. Nothing in it uploads. The language
bindings release afterwards from their own repositories; the last section
points to each one.

The tag commit is the source every package builds. A tag on a fork does
not publish: `.github/workflows/release.yml` runs only when
`github.repository == 'd-SEAMS/seams-core'`, with `GITHUB_TOKEN`. No other
token is added for this repository.

Build with `-j2`. Do not install into `/usr/lib`.

## 1. Tag

```sh
git checkout <release-commit>
git tag -a v2.14.0 -m "d-SEAMS 2.14.0"
git push https://github.com/d-SEAMS/seams-core.git v2.14.0
```

## 2. GitHub release

The tag push starts `.github/workflows/release.yml`, which cuts the
`## [2.14.0]` section out of `CHANGELOG.md` and runs `gh release create`
with it as the notes. Check the release page once the job finishes.

## 3. conda `libseams`

`recipe/recipe.yaml` builds `libseams` (`libyodaLib`, `libseams`,
`libseams_lang`, the `seams` CLI, the headers and the pkg-config and CMake
files) from the GitHub tag archive. Fill in the archive digest:

```sh
curl -L -o seams-core-2.14.0.tar.gz \
  https://github.com/d-SEAMS/seams-core/archive/refs/tags/v2.14.0.tar.gz
sha256sum seams-core-2.14.0.tar.gz
```

Write that digest into `recipe/recipe.yaml` (`source.sha256`). Then:

```sh
cat > variants.yaml <<'YAML'
c_compiler: [gcc]
c_compiler_version: ["14"]
cxx_compiler: [gxx]
cxx_compiler_version: ["14"]
c_stdlib: [sysroot]
c_stdlib_version: ["2.17"]
YAML
rattler-build build --recipe recipe -m variants.yaml -c conda-forge
```

`stdlib('c')` needs `c_stdlib` from a variant file. conda-forge's
pinning supplies it on the build farm; `variants.yaml` stands in for it
locally.

`rattler-build build` is the local proof. Copy `recipe/` into a
conda-forge/staged-recipes pull request (or the feedstock, once it
exists); opening that pull request is the publish step.

## 4. Yggdrasil `libseams_jll`

`packaging/yggdrasil/L/libseams/build_tarballs.jl` pins
`v2_14_0_tag_sha`. If `git rev-parse v2.14.0^{commit}` is not that
constant, replace it with the tag commit. Build locally:

```sh
cd packaging/yggdrasil/L/libseams
julia build_tarballs.jl --deploy=local
```

Copy `packaging/yggdrasil/L/libseams` into a Yggdrasil checkout and open
that pull request. When it merges, the registry bot registers
`libseams_jll` 2.14.0 in the Julia General registry, which Dseams.jl
needs before it can register.

## 5. Language bindings

Each binding moves its seams-core pin to `v2.14.0`, sets its version to
2.14.0, and releases by the steps its own repository documents:

| Repository | Package |
|------------|---------|
| [d-SEAMS/PydSEAMSlib](https://github.com/d-SEAMS/PydSEAMSlib) | `pydseamslib` on PyPI |
| [d-SEAMS/dseams-rs](https://github.com/d-SEAMS/dseams-rs) | crates `dseams-sys`, then `dseams` |
| [d-SEAMS/yodaStruct](https://github.com/d-SEAMS/yodaStruct) | LuaRocks `dseams` |
| [d-SEAMS/dseams-fortran](https://github.com/d-SEAMS/dseams-fortran) | fpm `dseams` |
| [d-SEAMS/dseams-r](https://github.com/d-SEAMS/dseams-r) | the R package `dseams` |
| [d-SEAMS/Dseams.jl](https://github.com/d-SEAMS/Dseams.jl) | `Dseams` in the General registry |

The first five can go in any order. Dseams.jl waits for `libseams_jll`.
The core Docs workflow installs `pydseamslib` from the PydSEAMSlib
`v2.14.0` tag, so the published core documentation needs that tag.
`docs/orgmode/contributing/repositories.org` has the pin locations and the
release order in one place.
