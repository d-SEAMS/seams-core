#!/bin/sh
# Copyright (c) d-SEAMS developers
# chemfiles 0.10.4 vendors fmt 6. fmt's fallback char8_t is an enum in
# namespace fmt. libc++ 19 leaves std::char_traits incomplete for that
# type, so basic_string_view<char8_t> does not instantiate. Count and
# compare the code units directly. The archive is repacked before
# chemfiles' external/CMakeLists.txt unpacks it.
set -eu
if [ ! -f external/fmt.tar.gz ]; then
  exit 0
fi
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
tar -xzf external/fmt.tar.gz -C "$tmp"
header=$tmp/fmt/include/fmt/core.h
if [ ! -f "$header" ]; then
  echo "chemfiles fmt archive has no include/fmt/core.h" >&2
  exit 1
fi
if grep -q 'fmt char length without char_traits' "$header"; then
  exit 0
fi
python3 - "$header" << 'PY'
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
text = path.read_text()
old_ctor = (
    "  basic_string_view(const Char* s)\n"
    "      : data_(s), size_(std::char_traits<Char>::length(s)) {}\n"
)
new_ctor = (
    "  // fmt char length without char_traits: libc++ 19 does not\n"
    "  // specialize std::char_traits for fmt's fallback char8_t.\n"
    "  basic_string_view(const Char* s) : data_(s), size_(0) {\n"
    "    while (s[size_] != 0)\n"
    "      ++size_;\n"
    "  }\n"
)
old_cmp = (
    "    int result = std::char_traits<Char>::compare(data_, other.data_, str_size);\n"
)
new_cmp = (
    "    int result = 0;\n"
    "    for (size_t i = 0; result == 0 && i < str_size; ++i) {\n"
    "      if (data_[i] < other.data_[i])\n"
    "        result = -1;\n"
    "      else if (other.data_[i] < data_[i])\n"
    "        result = 1;\n"
    "    }\n"
)
if old_ctor not in text or old_cmp not in text:
    raise SystemExit("fmt core.h does not match the 0.10.4 char_traits sites")
path.write_text(text.replace(old_ctor, new_ctor, 1).replace(old_cmp, new_cmp, 1))
PY
tar -czf external/fmt.tar.gz -C "$tmp" fmt
