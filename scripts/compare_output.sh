#!/usr/bin/env bash
# Run dfngen on every input in tests/data and compare the output with a baseline.
#
# Usage:
#   scripts/compare_output.sh record <baseline_dir>   # save the current output as baseline
#   scripts/compare_output.sh check  <baseline_dir>   # diff the current output against it
#
# DFN_output.txt is ignored because it contains timestamps.
set -euo pipefail

mode="${1:?usage: compare_output.sh record|check <baseline_dir>}"
baseline_dir="${2:?usage: compare_output.sh record|check <baseline_dir>}"
repo_root="$(cd "$(dirname "$0")/.." && pwd)"

cargo build --release --quiet --manifest-path "$repo_root/Cargo.toml"

case "$mode" in
record) output_dir="$baseline_dir" ;;
check) output_dir="$(mktemp -d)" ;;
*)
    echo "unknown mode: $mode" >&2
    exit 2
    ;;
esac

# Inputs refer to auxiliary files by relative path, so run from tests/data.
cd "$repo_root/tests/data"

status=0
for input_file in *.dat; do
    case_name="${input_file%.dat}"
    rm -rf "${output_dir:?}/$case_name"
    "$repo_root/target/release/dfngen" "$input_file" "$output_dir/$case_name" >/dev/null

    if [ "$mode" = check ]; then
        if diff -r -q -x DFN_output.txt "$baseline_dir/$case_name" "$output_dir/$case_name"; then
            echo "ok   $case_name"
        else
            echo "FAIL $case_name"
            status=1
        fi
    fi
done

exit "$status"
