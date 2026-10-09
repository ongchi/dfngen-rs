#!/usr/bin/env bash
# Run dfngen on every input in tests/data and compare the output with a baseline.
#
# Usage:
#   scripts/compare_output.sh record <baseline_dir>   # save the current output as baseline
#   scripts/compare_output.sh check  <baseline_dir>   # the output must be identical to it
#   scripts/compare_output.sh check  <baseline_dir> <tolerance>
#                                                     # decimal numbers can differ by the tolerance
#
# The console log is compared with the timestamps removed.
set -euo pipefail

mode="${1:?usage: compare_output.sh record|check <baseline_dir>}"
baseline_dir="${2:?usage: compare_output.sh record|check <baseline_dir>}"
tolerance="${3:-}"
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
mkdir -p "$output_dir"

cd "$repo_root/tests/data"

status=0
for input_file in *.toml; do
    case_name="${input_file%.toml}"
    rm -rf "${output_dir:?}/$case_name"
    # Keep the console log too, without the leading timestamp of each line
    # and without the run-specific output path.
    RUST_LOG=info "$repo_root/target/release/dfngen" "$input_file" "$output_dir/$case_name" |
        sed -E -e 's/^[^ ]+ //' -e "s#$output_dir#<output>#g" >"$output_dir/$case_name.log"
    mv "$output_dir/$case_name.log" "$output_dir/$case_name/console.log"

    if [ "$mode" = check ]; then
        if diff -r -q"$baseline_dir/$case_name" "$output_dir/$case_name" >/dev/null; then
            echo "ok   $case_name (identical)"
        elif [ -n "$tolerance" ] && python3 "$repo_root/scripts/compare_output.py" \
            "$baseline_dir/$case_name" "$output_dir/$case_name" "$tolerance"; then
            echo "ok   $case_name (within tolerance $tolerance)"
        else
            echo "FAIL $case_name (output kept in $output_dir/$case_name)"
            status=1
        fi
    fi
done

exit "$status"
