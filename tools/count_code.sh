#!/usr/bin/env bash
# Counts kentron's own code with cloc: what was written for the project, not generated.
#
# Left out:
#   - Rust files marked `@generated` on their first line (tests written by tools/draw2rs.py);
#   - YAML cases under crates/occt-tests/cases/occt/ (written by tools/drawyaml.py and
#     tools/gtest2yaml.py);
#   - the not-implemented command stubs tools/draw2rs.py appends to
#     crates/occt-tests/src/draw/commands.rs (the implemented commands there count);
#   - untracked files (scratch examples), documentation and JSON data.
# Counted: the library with its unit tests, the hand-ported test suites, the test harness,
# pykentron and the tools: one table, a row per part and the total, code lines by language.
#
# Usage: tools/count_code.sh [--by-file]   (--by-file: also every file's counts)

set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
command -v cloc >/dev/null || { echo "cloc is not installed" >&2; exit 1; }

by_file=()
[[ "${1:-}" == "--by-file" ]] && by_file=(--by-file)

languages="Rust,Python,YAML,Bourne Shell,Bourne Again Shell,TOML"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# commands.rs without the generated stubs: a run of `///` lines, then
#   pub fn name(&mut self, args: &[&str]) -> DrawResult {
#       self.not_implemented("name", args)
#   }
commands=crates/occt-tests/src/draw/commands.rs
stripped="$work/commands.rs"
awk '
    { line[NR] = $0 }
    END {
        i = 1
        while (i <= NR) {
            j = i
            while (j <= NR && line[j] ~ /^[ \t]*\/\/\//) j++
            if (j + 2 <= NR && line[j] ~ /pub fn [A-Za-z0-9_]+\(&mut self, args: &\[&str\]\) -> DrawResult \{/ \
                && line[j + 1] ~ /self\.not_implemented\(/ && line[j + 2] ~ /^[ \t]*\}[ \t]*$/) {
                i = j + 3
                continue
            }
            print line[i]
            i++
        }
    }' "$commands" > "$stripped"

# The authored files among those git tracks under the given paths, one per line.
authored() {
    git ls-files -- "$@" | while IFS= read -r f; do
        case "$f" in
            crates/occt-tests/cases/occt/*) continue ;;
            "$commands") echo "$stripped"; continue ;;
        esac
        if [[ "$f" == *.rs ]] && head -n 1 "$f" | grep -q '@generated'; then
            continue
        fi
        echo "$f"
    done
}

# One row: a part's files, code lines by language, comment lines and code lines in all.
row() {
    local title=$1 list=$2
    cloc --quiet --csv --include-lang="$languages" --list-file="$list" 2>/dev/null | awk -F, -v title="$title" '
        NR == 1 { next }
        $2 == "SUM" { files = $1; comment = $4; code = $5; next }
        $2 == "Rust" { rust = $5; next }
        $2 == "Python" { python = $5; next }
        $2 == "YAML" { yaml = $5; next }
        { other += $5 }
        END { printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", title, files, rust, python, yaml, other, comment, code }'
}

authored crates/kentron > "$work/library"
authored crates/occt-tests/src crates/occt-tests/Cargo.toml > "$work/harness"
authored crates/occt-tests/tests crates/occt-tests/cases > "$work/tests"
authored crates/pykentron > "$work/pykentron"
authored tools Cargo.toml > "$work/tools"
cat "$work"/library "$work"/harness "$work"/tests "$work"/pykentron "$work"/tools > "$work/all"

{
    row "kentron library (with its unit tests)" "$work/library"
    row "OCCT test harness" "$work/harness"
    row "tests written for kentron" "$work/tests"
    row "pykentron (bindings, CadQuery layer)" "$work/pykentron"
    row "tools and workspace" "$work/tools"
    row "total" "$work/all"
} | awk -F'\t' '
    function line() { printf "%s\n", "----------------------------------------+-------+---------+---------+--------+-------+----------+---------" }
    function fmt(n) { s = sprintf("%d", n); out = ""; while (length(s) > 3) { out = "," substr(s, length(s) - 2) out; s = substr(s, 1, length(s) - 3) }; return s out }
    BEGIN {
        printf "%-40s| %5s | %7s | %7s | %6s | %5s | %8s | %7s\n", "Part", "Files", "Rust", "Python", "YAML", "Other", "Comments", "Code"
        line()
    }
    $1 == "total" { line() }
    { printf "%-40s| %5s | %7s | %7s | %6s | %5s | %8s | %7s\n", $1, fmt($2), fmt($3), fmt($4), fmt($5), fmt($6), fmt($7), fmt($8) }'
echo "(lines of code by language; Other: shell scripts and TOML; blank lines not counted)"

if ((${#by_file[@]})); then
    echo
    cloc --quiet --hide-rate --include-lang="$languages" --by-file --list-file="$work/all" | sed '/^$/d'
fi

generated=$(git ls-files -- '*.rs' | while IFS= read -r f; do if head -n 1 "$f" | grep -q '@generated'; then echo "$f"; fi; done | wc -l)
yaml=$(git ls-files -- 'crates/occt-tests/cases/occt/*' | wc -l)
stubs=$(( $(grep -c 'pub fn' "$commands") - $(grep -c 'pub fn' "$stripped") ))
echo "Left out: $generated generated Rust files, $yaml generated YAML files, $stubs generated command stubs."
