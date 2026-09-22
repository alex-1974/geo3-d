#!/usr/bin/env python3
"""
Print D import paths resolved by DUB for the current package graph.

One absolute path is written per line. Consumers can translate the
result to compiler-specific `-I...` arguments without knowing dependency
locations themselves.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def die(message: str) -> None:
    raise SystemExit(
        f"error: {message}"
    )


def find_repo_root() -> Path:
    return (
        Path(__file__)
        .resolve()
        .parent
        .parent
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument(
        "--compiler",
        default="dmd",
        help=(
            "D compiler passed to "
            "`dub describe` "
            "(default: dmd)"
        ),
    )

    args = parser.parse_args()

    root = find_repo_root()

    result = subprocess.run(
        [
            "dub",
            "describe",
            "--data=import-paths",
            "--data-list",
            f"--compiler={args.compiler}",
        ],
        cwd=root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    if result.returncode != 0:
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)

        die("`dub describe` failed")

    resolved: list[Path] = []
    seen: set[Path] = set()

    for line in result.stdout.splitlines():
        value = line.strip()

        if not value:
            continue

        path = Path(value).resolve()

        if not path.is_dir():
            die(
                "resolved import path "
                "is not a directory: "
                + str(path)
            )

        if path in seen:
            continue

        seen.add(path)
        resolved.append(path)

    if not resolved:
        die("DUB resolved no import paths")

    for path in resolved:
        print(path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
