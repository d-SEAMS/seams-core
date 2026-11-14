#!/bin/bash
set -euxo pipefail

# system_libseams builds only the extension, linked against the
# libseams output in the host prefix.
"${PYTHON}" -m pip install . -vv --no-deps --no-build-isolation \
  -Csetup-args=-Dsystem_libseams=true \
  -Csetup-args=--wrap-mode=nodownload \
  -Cbuild-dir=_build_py
