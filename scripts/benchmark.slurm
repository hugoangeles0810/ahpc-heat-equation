#!/usr/bin/env bash
#SBATCH --job-name=heat-bench
#SBATCH --output=results/slurm-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=00:30:00
# ----------------------------------------------------------------------------
# Strong-scaling sweep for the MPI heat solver.
#
# Submit with:   sbatch scripts/benchmark.slurm
#
# For each process count P in PROCS the job runs the solver and appends the
# parsed timing to results/benchmark.csv via collect_results.sh. The raw
# stdout of every run is kept under results/raw/ for traceability.
#
# NOTE: the current solver hardcodes the grid size (imax=kmax=80). To also
# sweep the problem size n (objective b), first parameterize the grid from
# argv in src/heat-mpi-big.c, then add an N loop around the P loop below and
# pass $N to the binary.
# ----------------------------------------------------------------------------
set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

PROCS=(1 2 4 8 16)
BIN="bin/heat-mpi-big"
CSV="results/benchmark.csv"
RAW="results/raw"

module load openmpi 2>/dev/null || true   # adjust to your cluster's module name

make mpi
mkdir -p "$RAW"

# CSV header (write once)
if [[ ! -f "$CSV" ]]; then
  echo "nprocs,idim,kdim,iterations,wall_time_s,comm_time_s,criterion_time_s" > "$CSV"
fi

for P in "${PROCS[@]}"; do
  OUT="$RAW/run_p${P}.out"
  echo ">> P=$P"
  mpirun -np "$P" "$BIN" | tee "$OUT"
  scripts/collect_results.sh "$OUT" >> "$CSV"
done

echo ">> wrote $CSV"
