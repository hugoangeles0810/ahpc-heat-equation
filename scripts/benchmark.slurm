#!/usr/bin/env bash
#SBATCH --job-name=heat-bench
#SBATCH --output=results/slurm-%j.out
#SBATCH --partition=standard
#SBATCH --nodelist=n003,n005
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks-per-socket=4
#SBATCH --hint=nomultithread
#SBATCH --time=01:00:00
# Plan B: two-node strong scaling on n003+n005 (homogeneous Xeon Gold 6130),
# 8 physical cores per node (4 per socket). Total 16 cores = 32 logical CPUs,
# exactly the account cap (cpu=32). --hint=nomultithread gives whole cores;
# --ntasks-per-socket=4 forces a balanced 4-cores-per-NUMA-domain layout so
# the ppr placement below has both sockets to fill. No --exclusive (would
# grab all 64 logical CPUs per node -> rejected with QOSMaxCpuPerUserLimit).
# vs exp01 (16 ranks on ONE node): here P=16 is 8+8 across two nodes, so the
# 8->16 step now crosses the network and isolates inter-node MPI cost.
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
PROCS=(1 2 4 8 16)            # MPI procs (1 rank/core; 1..8 on n003, 16 = 8+8 across n003+n005)
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
    if mpirun --map-by ppr:4:socket --bind-to core -np "$P" "$BIN" "$N" "$ITMAX" 0 | tee "$OUT"; then
      scripts/collect_results.sh "$OUT" >> "$CSV"
    else
      echo "!! run N=$N P=$P FAILED, skipping" >&2
    fi
  done
done

echo ">> wrote $CSV"
