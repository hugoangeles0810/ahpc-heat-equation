#!/usr/bin/env bash
# Parse the solver's summary line into one CSV row.
#
#   scripts/collect_results.sh run.out      # read from a file
#   mpirun ... | scripts/collect_results.sh # read from stdin
#
# The MPI solver prints a machine-readable summary line, e.g.
#   !     4 =  2 x  2  12345      0.1234       1       0.0123     1     10     0.0045
# i.e.  !<nprocs> =<idim> x<kdim> <iters> <wall> 1 <comm> 1 <stride> <criterion>
#
# Output columns (no header — benchmark.slurm writes the header once):
#   nprocs,idim,kdim,iterations,wall_time_s,comm_time_s,criterion_time_s
set -euo pipefail

awk '
  /^!/ {
    line = $0
    gsub(/[!=x]/, " ", line)
    n = split(line, f)
    if (f[1] ~ /^[0-9]+$/)           # skip the two header lines (non-numeric first field)
      print f[1] "," f[2] "," f[3] "," f[4] "," f[5] "," f[7] "," f[10]
  }
' "${1:-/dev/stdin}"
