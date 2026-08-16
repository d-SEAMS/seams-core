# Contributing

We love pull requests from everyone. By participating in this project, you
agree to abide by the contributor covenant suggested [code of conduct].

[code of conduct]: https://github.com/d-SEAMS/seams-core/blob/master/CODE_OF_CONDUCT.md

Do make sure to run tests and generally **BE PREPARED** to have your code vetted and
checked. Do **NOT** submit code you would not be prepared to defend or maintain,
unless you are fixing a bug.

Push to your fork and [submit a pull request][pr].

[pr]: https://github.com/d-SEAMS/seams-core/compare

At this point you're waiting on us. We like to at least comment on pull requests
within four business days (and, typically, three business day). We may suggest
some changes or improvements or alternatives.

Some things that will increase the chance that your pull request is accepted:

- Write tests (Catch2; `pixi run test` or `meson test -C bbdir`).
- Follow Conventional Commits below.

## Commit style

Subjects use Conventional Commits:

```
feat(cli): add --print-config for the twelve-factor table
fix(neighbours): take the LAMMPS dump 3x3 for cutoff lists
docs(book): document vesin versus linkcell
```

Prefixes that land here: `feat`, `fix`, `docs`, `chore`, `refactor`,
`test`. Scope is optional (`cli`, `gpu`, `book`). The 2020
`fileName: ...` sandwich sample is not the house style.

Name every co-author with a trailer when the work is shared:

```
Co-authored-by: Name <email>
```

Build with pixi or the Nix flake. The engine binary is `seams`
(`./bbdir/seams` after `pixi run build`). Python is PydSEAMSlib;
Lua is yodaStruct. Do not add a `yodaStruct -c` / `conf.yaml`
driver here.