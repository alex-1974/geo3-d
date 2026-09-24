#!/usr/bin/env python3

from pathlib import Path


ROOT = Path(__file__).resolve().parent

SOURCE = (
    ROOT
    / ".."
    / ".."
    / "source"
    / "geo3"
    / "internal"
    / "orientation_filter.d"
).resolve()

DESTINATION = (
    ROOT
    / "source"
    / "production_filter_local.d"
)


SOURCE_MODULE = "module geo3.internal.orientation_filter;"
LOCAL_MODULE = "module production_filter_local;"


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")

    count = text.count(SOURCE_MODULE)

    if count != 1:
        raise RuntimeError(
            "expected exactly one production orientation-filter "
            f"module declaration, found {count}"
        )

    generated = text.replace(
        SOURCE_MODULE,
        LOCAL_MODULE,
        1,
    )

    DESTINATION.write_text(
        generated,
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
