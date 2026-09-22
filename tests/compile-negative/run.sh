#!/usr/bin/env bash
set -euo pipefail

DC="${1:-${DC:-dmd}}"

ROOT="$(
    cd "$(dirname "${BASH_SOURCE[0]}")/../.."
    pwd
)"

FIXTURES="$ROOT/tests/compile-negative"
TMP="$(mktemp -d)"

trap 'rm -rf "$TMP"' EXIT

cd "$ROOT"

mapfile -t IMPORT_PATHS < <(
    PYTHONDONTWRITEBYTECODE=1 \
        python3 \
        "$ROOT/tools/dub-import-paths.py" \
        --compiler "$DC"
)

IMPORT_FLAGS=()

for import_path in "${IMPORT_PATHS[@]}"; do
    IMPORT_FLAGS+=(
        "-I$import_path"
    )
done

failures=0


compile_positive()
{
    family="$1"
    source="$FIXTURES/${family}_positive.d"
    object="$TMP/${family}_positive.o"
    log="$TMP/${family}_positive.log"

    if "$DC" \
        -preview=dip1000 \
        "${IMPORT_FLAGS[@]}" \
        -c "$source" \
        "-of=$object" \
        >"$log" 2>&1
    then
        printf 'PASS positive: %s\n' "$family"
    else
        printf 'FAIL positive: %s\n' "$family"
        cat "$log"
        failures=$((failures + 1))
    fi
}


compile_negative()
{
    family="$1"
    source="$FIXTURES/${family}_negative.d"
    object="$TMP/${family}_negative.o"
    log="$TMP/${family}_negative.log"

    if "$DC" \
        -preview=dip1000 \
        "${IMPORT_FLAGS[@]}" \
        -c "$source" \
        "-of=$object" \
        >"$log" 2>&1
    then
        printf 'FAIL negative: %s unexpectedly compiled\n' "$family"
        failures=$((failures + 1))
    else
        printf 'PASS negative: %s rejected\n' "$family"
    fi
}


echo "compiler: $DC"
"$DC" --version | sed -n '1,3p'
echo

for family in polyline ring
do
    compile_positive "$family"
    compile_negative "$family"
done

echo
echo "lifetime compile-negative failures: $failures"

if (( failures != 0 )); then
    exit 1
fi
