#!/usr/bin/env python3
"""Compare two dfngen output folders with a numeric tolerance.

Usage: compare_output.py <baseline_dir> <output_dir> <tolerance>

The folders must contain the same files with the same number of lines and the
same number of values on each line. Integers and text must be identical.
Decimal numbers must agree within the tolerance (relative, or absolute for
values smaller than 1).

Exit status: 0 if the folders agree, 1 if not.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

IGNORED_FILES = {"DFN_output.txt"}
TOKEN_SEPARATORS = re.compile(r"[\s,{}()\[\]]+")
MAX_REPORTED_DIFFERENCES = 5


def list_files(root: Path) -> set[Path]:
    return {
        path.relative_to(root)
        for path in root.rglob("*")
        if path.is_file() and path.name not in IGNORED_FILES
    }


def parse_decimal(token: str) -> float | None:
    """Return the value of a decimal number, or None for integers and text."""
    if token.lstrip("+-").isdigit():
        return None
    try:
        return float(token)
    except ValueError:
        return None


def compare_tokens(expected: str, actual: str, tolerance: float) -> float | None:
    """Return the difference of two tokens, or None if they do not agree."""
    if expected == actual:
        return 0.0

    expected_value = parse_decimal(expected)
    actual_value = parse_decimal(actual)
    if expected_value is None or actual_value is None:
        return None

    difference = abs(expected_value - actual_value)
    scale = max(1.0, abs(expected_value), abs(actual_value))
    return difference / scale if difference <= tolerance * scale else None


def compare_files(expected_path: Path, actual_path: Path, tolerance: float):
    """Return (largest difference, list of messages for values out of tolerance)."""
    expected_lines = expected_path.read_text().splitlines()
    actual_lines = actual_path.read_text().splitlines()

    if len(expected_lines) != len(actual_lines):
        return 0.0, [f"{len(expected_lines)} lines expected, {len(actual_lines)} found"]

    largest_difference = 0.0
    messages = []

    for line_number, (expected_line, actual_line) in enumerate(
        zip(expected_lines, actual_lines), start=1
    ):
        expected_tokens = TOKEN_SEPARATORS.split(expected_line.strip())
        actual_tokens = TOKEN_SEPARATORS.split(actual_line.strip())

        if len(expected_tokens) != len(actual_tokens):
            messages.append(f"line {line_number}: different number of values")
            continue

        for expected, actual in zip(expected_tokens, actual_tokens):
            difference = compare_tokens(expected, actual, tolerance)
            if difference is None:
                messages.append(f"line {line_number}: expected {expected}, found {actual}")
            else:
                largest_difference = max(largest_difference, difference)

    return largest_difference, messages


def main() -> int:
    if len(sys.argv) != 4:
        print(__doc__, file=sys.stderr)
        return 2

    baseline_dir, output_dir = Path(sys.argv[1]), Path(sys.argv[2])
    tolerance = float(sys.argv[3])

    baseline_files = list_files(baseline_dir)
    output_files = list_files(output_dir)
    status = 0

    for path in sorted(baseline_files ^ output_files):
        owner = "baseline" if path in baseline_files else "output"
        print(f"  {path}: only in {owner}")
        status = 1

    largest_difference = 0.0
    changed_files = 0

    for path in sorted(baseline_files & output_files):
        difference, messages = compare_files(baseline_dir / path, output_dir / path, tolerance)
        largest_difference = max(largest_difference, difference)
        changed_files += difference > 0.0 or bool(messages)

        if messages:
            status = 1
            print(f"  {path}: {len(messages)} differences out of tolerance")
            for message in messages[:MAX_REPORTED_DIFFERENCES]:
                print(f"    {message}")

    print(
        f"  {changed_files} of {len(baseline_files & output_files)} files differ,"
        f" largest difference in tolerance: {largest_difference:.3e}"
    )
    return status


if __name__ == "__main__":
    sys.exit(main())
