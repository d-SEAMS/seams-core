"""The cutoff-list page records this tip's CSR session."""

import sys
from pathlib import Path

PAGE = Path("docs/orgmode/howto/neighbour-lists.org")


def main() -> int:
    text = PAGE.read_text(encoding="utf-8")
    required = (
        "one session",
        "8 cores",
        "interquartile range (IQR)",
        "5th of the 8 sorted samples",
        "No host name",
        "4.554",
        "6.879",
        "5.136",
        "0.322",
        "1.388",
        "0.206",
    )
    missing = [item for item in required if item not in text]
    stale = [item for item in ("5.128", "6.093", "4.882") if item in text]
    if missing or stale:
        sys.stderr.write(
            "neighbour-list record missing "
            + ", ".join(missing)
            + " stale "
            + ", ".join(stale)
            + "\n"
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
