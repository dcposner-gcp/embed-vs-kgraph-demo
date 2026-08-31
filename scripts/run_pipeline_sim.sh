#!/usr/bin/env bash
# scripts/run_pipeline_sim.sh -- run the UNCHANGED targets pipeline against a
# simulated data directory without touching data/ or the project's _targets/
# store. Builds scratch/sim_run/<name>/ with symlinks to _targets.R and R/,
# links data -> <sim dir>, and runs tar_make() there with its own store,
# scratch/ and output/. Host R with the project's renv library (falls back to
# whatever `Rscript` resolves, e.g. inside the container).
#
#   scripts/run_pipeline_sim.sh [SIM_DIR=data_sim/phecode] [NAME=basename(SIM_DIR)]
#
# Outputs: scratch/sim_run/<name>/output/dx/  (+ _targets/ store, scratch/ QC files)
set -euo pipefail

SIM_DIR="${1:-data_sim/phecode}"
NAME="${2:-$(basename "$SIM_DIR")}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SIM_ABS="$(cd "$SIM_DIR" && pwd)"
RUN="$ROOT/scratch/sim_run/$NAME"

[[ -f "$SIM_ABS/validation/dict-disease-disease-wiki-bin-eval.Rdata" ]] || {
  echo "ERROR: $SIM_DIR does not look like a simulated data tree (run scripts/simulate_data.R first)" >&2; exit 1; }

mkdir -p "$RUN"
ln -sfn "$ROOT/_targets.R" "$RUN/_targets.R"
ln -sfn "$ROOT/R"          "$RUN/R"
ln -sfn "$SIM_ABS"         "$RUN/data"

LIB="$(ls -d "$ROOT"/renv/library/*/R-*/* 2>/dev/null | head -n 1 || true)"
cd "$RUN"
echo "run dir : $RUN"
echo "data    : $SIM_ABS"
if [[ -n "$LIB" && -d "$LIB" ]]; then
  echo "R lib   : $LIB"
  R_PROFILE_USER=/dev/null R_LIBS="$LIB${R_LIBS:+:$R_LIBS}" \
    Rscript --vanilla -e 'targets::tar_make(reporter = "balanced")'
else
  echo "R lib   : (default library path; no renv library found under $ROOT/renv/library)"
  Rscript -e 'targets::tar_make(reporter = "balanced")'
fi
echo
echo "outputs : $RUN/output/dx/"
ls -1 "$RUN/output/dx/" 2>/dev/null || true
