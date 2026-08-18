#!/usr/bin/env bash
# End-to-end comparison against the reference R implementation.
# Run from this directory on a machine that has R.
#
#   ./run_comparison.sh [replicates]     (default 5000)
set -euo pipefail

REPLICATES="${1:-5000}"

# --- locate the interpreters ------------------------------------------------
# Neither name is reliably on PATH, so prefer an explicit override, then PATH,
# then the usual Windows install location. Python candidates are probed by
# RUNNING them rather than by resolving the name: on Windows "python3" normally
# resolves to a Microsoft Store stub that prints an install message and does
# nothing, which command -v happily reports as success.
PY="${PYTHON:-}"
if [ -z "$PY" ]; then
    for candidate in python3 python py; do
        if command -v "$candidate" >/dev/null 2>&1 &&
           "$candidate" -c "import numpy" >/dev/null 2>&1; then
            PY="$candidate"
            break
        fi
    done
fi
[ -n "$PY" ] || {
    echo "No working Python with numpy found. Set PYTHON=/path/to/python." >&2
    exit 1
}

RSCRIPT="${RSCRIPT:-}"
if [ -z "$RSCRIPT" ] && command -v Rscript >/dev/null 2>&1; then
    RSCRIPT="Rscript"
fi
if [ -z "$RSCRIPT" ]; then
    # Windows: take the highest version number installed under Program Files.
    RSCRIPT="$(ls -d "/c/Program Files/R/R-"*/bin/Rscript.exe 2>/dev/null | sort -V | tail -1 || true)"
fi
[ -n "$RSCRIPT" ] || {
    echo "Rscript not found. Install R, then inside R:  install.packages('BalancedSampling')" >&2
    echo "If R is installed somewhere unusual, set RSCRIPT=/path/to/Rscript." >&2
    exit 1
}

echo "python:  $("$PY" --version 2>&1)"
echo "R:       $("$RSCRIPT" -e 'cat(R.version.string)' 2>/dev/null | tail -1)"

"$RSCRIPT" -e 'if (!requireNamespace("BalancedSampling", quietly=TRUE)) { cat("Installing BalancedSampling...\n"); install.packages("BalancedSampling", repos="https://cloud.r-project.org") }'

echo "== building the shared frame"
"$PY" frames.py

echo "== deterministic checks: same inputs, so the outputs must match exactly"
"$PY" diagnostic.py
"$RSCRIPT" diagnostic.R

echo "== python draws"
"$PY" run_python.py --replicates "$REPLICATES" --method lpm2
"$PY" run_python.py --replicates "$REPLICATES" --method lpm1

echo "== R draws"
"$RSCRIPT" run_r.R "$REPLICATES" lpm2
"$RSCRIPT" run_r.R "$REPLICATES" lpm1

echo "== comparison"
"$PY" compare.py
