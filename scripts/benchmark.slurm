#!/usr/bin/env bash
#SBATCH --job-name=heat-bench
#SBATCH --output=results/slurm-%j.out
#SBATCH --partition=standard
#SBATCH --nodelist=n003,n004,n005
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks=48
#SBATCH --exclusive
#SBATCH --time=01:00:00
# ----------------------------------------------------------------------------
# Scaling sweep for the MPI heat solver.
#
# Submit with:   sbatch scripts/benchmark.slurm
#
# The solver takes its parameters from argv:
#   bin/heat-mpi-big <grid_size> <itmax> <print_result>
#
# This job sweeps two axes:
#   - problem size N  (NS array)   -> grid is imax=kmax=N
#   - process count P (PROCS array)
# giving both strong scaling (fixed N, varying P) and problem-size scaling.
#
# For each (N, P) pair the solver runs and its summary line is parsed by
# collect_results.sh into one CSV row of results/benchmark.csv. The raw stdout
# of every run is kept under results/raw/ for traceability.
# ----------------------------------------------------------------------------
set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

NS=(128 256 512 1024)         # grid sizes (imax=kmax=N)
PROCS=(1 2 4 8 16 32 48)      # MPI process counts (1..16 = 1 node; 32 = 2 nodes; 48 = 3 nodes)
ITMAX=20000                   # iteration cap (solver stops earlier on convergence)
BIN="bin/heat-mpi-big"
CSV="results/benchmark.csv"
RAW="results/raw"

module load gnu12/12.4.0
module load openmpi4/4.1.6

make mpi
mkdir -p "$RAW"

# CSV header (write once)
if [[ ! -f "$CSV" ]]; then
  echo "grid_size,nprocs,idim,kdim,iterations,wall_time_s,comm_time_s,criterion_time_s" > "$CSV"
fi

for N in "${NS[@]}"; do
  for P in "${PROCS[@]}"; do
    OUT="$RAW/run_n${N}_p${P}.out"
    echo ">> N=$N P=$P"
    if mpirun --map-by ppr:8:socket --bind-to core -np "$P" "$BIN" "$N" "$ITMAX" 0 | tee "$OUT"; then
      scripts/collect_results.sh "$OUT" >> "$CSV"
    else
      echo "!! run N=$N P=$P FAILED, skipping" >&2
    fi
  done
done

echo ">> wrote $CSV"
