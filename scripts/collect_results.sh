#!/usr/bin/env bash
# Extract the solver's RESULT row into CSV.
#
#   scripts/collect_results.sh run.out      # read from a file
#   mpirun ... | scripts/collect_results.sh # read from stdin
#
# Both solvers (serial heat-big and MPI heat-mpi-big) print one identical,
# self-describing summary line; the serial run reports nprocs=idim=kdim=1 and
# comm=criterion=0:
#   RESULT,<grid_size>,<nprocs>,<idim>,<kdim>,<iterations>,<wall_s>,<comm_s>,<crit_s>
#
# Collecting it is just "keep the RESULT lines, drop the RESULT, tag".
#
# Output columns (no header — benchmark.slurm writes the header once):
#   grid_size,nprocs,idim,kdim,iterations,wall_time_s,comm_time_s,criterion_time_s
set -euo pipefail

awk '/^RESULT,/ { print substr($0, 8) }' "${1:-/dev/stdin}"
