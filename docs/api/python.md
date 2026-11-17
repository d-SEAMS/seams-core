# Python

The Python binding is the `pydseamslib` distribution, with `import dseams`.
This release's distribution is being published to PyPI from its own repository,
[d-SEAMS/PydSEAMSlib](https://github.com/d-SEAMS/PydSEAMSlib), with its API reference, tutorials,
tests and release notes. It links the engine through the C API
(`libseams`, `include/seams.h`) and the C ABI `libseams_lang`
(`bindings/cxx/seams_lang.h`) that this repository builds and installs.

The binding documentation is at <https://python.dseams.info/>, and the Python walkthrough is [on the video page](https://docs.dseams.info/videos/embeds/python.html).

```bash
pip install 'pydseamslib==@PROJECT_VERSION@'
```

Until that upload is visible on the selected package index, install the
binding from source against this checkout:

```bash
sh scripts/install_pydseamslib.sh
```

`docs/orgmode/contributing/repositories.org` describes how the repositories fit together.
