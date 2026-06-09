#!/usr/bin/env bash
# Quick local smoke test of the MPI solver on your laptop.
#
#   scripts/run_local.sh            # run with 4 ranks (default)
#   scripts/run_local.sh 2          # run with 2 ranks
#
# Builds first if the binary is missing. Use --oversubscribe so you can
# request more ranks than physical cores while developing.
set -euo pipefail

cd "$(dirname "$0")/.."

NPROCS="${1:-4}"
BIN="bin/heat-mpi-big"

if [[ ! -x "$BIN" ]]; then
  echo ">> building ($BIN not found)"
  make mpi
fi

echo ">> mpirun -np $NPROCS $BIN"
mpirun --oversubscribe -np "$NPROCS" "$BIN"
