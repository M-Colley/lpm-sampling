#!/usr/bin/env bash
# End-to-end comparison against the reference R implementation.
# Run from this directory on a machine that has R.
set -euo pipefail

REPLICATES="${1:-5000}"

command -v Rscript >/dev/null || {
    echo "Rscript not found. Install R, then inside R:  install.packages('BalancedSampling')" >&2
    exit 1
}

Rscript -e 'if (!requireNamespace("BalancedSampling", quietly=TRUE)) { cat("Installing BalancedSampling...\n"); install.packages("BalancedSampling", repos="https://cloud.r-project.org") }'

echo "== building the shared frame"
PYTHONPATH=../src python3 frames.py

echo "== python draws"
python3 run_python.py --replicates "$REPLICATES" --method lpm2
python3 run_python.py --replicates "$REPLICATES" --method lpm1

echo "== R draws"
Rscript run_r.R "$REPLICATES" lpm2
Rscript run_r.R "$REPLICATES" lpm1

echo "== comparison"
python3 compare.py
