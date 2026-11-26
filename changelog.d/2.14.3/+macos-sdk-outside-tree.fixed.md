The macOS CI SDK archive is downloaded under `RUNNER_TEMP` instead of the source tree, so the checkout does not gain an untracked `MacOSX11.0.sdk.tar.xz` before the source-tree tests.
