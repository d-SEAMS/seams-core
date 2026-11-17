An `unsigned-char` CI job builds every target on x86_64 Linux with `-funsigned-char` and runs the suite on it, so code that reads a plain `char` as signed fails without an aarch64 runner.
