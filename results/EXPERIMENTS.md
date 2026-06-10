# Experiments Log

Tracking of MPI heat-solver benchmark runs on the **Khipu** cluster (UTEC).
Each experiment lives in its own folder with the raw `benchmark.csv` and the four
generated figures (`time` / `speedup` / `efficiency` / `iterations`).

Folder naming: `expNN-<nodes>node-<cores-per-node>core`.

## Fixed setup (constant across all experiments unless noted)

- **Account:** `postgrado` — hard cap `cpu=32` *logical* CPUs per user.
- **Nodes:** homogeneous Intel Xeon Gold 6130 (`n003`, `n005`), 2 sockets × 16 cores,
  hyperthreading **on** (32 cores = 64 logical CPUs/node).
- **CPU policy:** 1 MPI rank per **physical** core (`--hint=nomultithread`), **no**
  `--exclusive` (would allocate all 64 logical CPUs/node → `QOSMaxCpuPerUserLimit`).
  Usable budget is therefore **16 physical cores** total.
- **Solver:** `bin/heat-mpi-big <N> <ITMAX> 0`, `ITMAX=20000`, GCC 12.4.0 + OpenMPI 4.1.6.
- **Job script:** [`scripts/benchmark.slurm`](../scripts/benchmark.slurm) is re-tuned per
  experiment; the per-experiment SBATCH/placement config is recorded below.

## Index

| Exp | Folder | Date | Nodes × cores | NS (grid sizes) | PROCS | One-line takeaway |
|-----|--------|------|---------------|-----------------|-------|-------------------|
| 01 | [exp01-1node-16core](exp01-1node-16core/) | 2026-06-10 | 1 × 16 (n003) | 128, 256, 512, 1024 | 1, 2, 4, 8, 16 | Baseline single-node strong scaling |
| 02 | [exp02-2node-8core](exp02-2node-8core/) | 2026-06-10 | 2 × 8 (n003+n005) | 128, 256, 512, 1024 | 1, 2, 4, 8, 16 | P=16 = 8+8 crosses the network; **beats** 1-node |

## Headline comparison (P=16, 16 ranks total)

Same 16 ranks, packed on one node (exp01) vs split 8+8 across two nodes (exp02):

| N | exp01 wall / eff / comm | exp02 wall / eff / comm |
|---|-------------------------|-------------------------|
| 512  | 0.973 s / 0.77 / 0.205 s | **0.890 s / 0.84 / 0.195 s** |
| 1024 | 3.472 s / 0.84 / 0.459 s | **3.205 s / 0.94 / 0.351 s** |

**Finding:** spreading ranks across two nodes is *faster* than consolidating them — the
solver is memory-bandwidth bound, so 2 nodes' aggregate memory bandwidth outweighs the
inter-node halo cost (and `comm_time` is even lower with only 8 ranks/node). Counter to
the naive "avoid the network" intuition.

## Per-experiment configuration

### exp01 — single node, 16 cores
- **SBATCH:** `--nodelist=n003 --nodes=1 --ntasks=16 --hint=nomultithread`
- **Placement:** `mpirun --map-by socket --bind-to core` (cyclic across the 2 sockets)
- **Result:** [exp01-1node-16core/](exp01-1node-16core/) — strong scaling per grid size;
  small grids (n=128) saturate early (S≈4 at p≥8) as halo comm dominates; large grids
  scale near-ideal (n=1024: S=13.4 at p=16).

### exp02 — two nodes, 8 cores each
- **SBATCH:** `--nodelist=n003,n005 --nodes=2 --ntasks=16 --ntasks-per-node=8 --ntasks-per-socket=4 --hint=nomultithread`
- **Placement:** `mpirun --map-by ppr:4:socket --bind-to core` (block-fill n003 for P≤8,
  then 8+8 across n003+n005 at P=16)
- **Result:** [exp02-2node-8core/](exp02-2node-8core/) — see headline comparison above.
- **Note:** intra-node placement is *block* (`ppr`) here vs *cyclic* (`--map-by socket`)
  in exp01, so only the P=16 (all-cores) point is a strictly apples-to-apples comparison.

## Planned / next experiments

- Vary **NS** (problem size) — e.g. larger grids (2048, 4096) to push further into the
  memory-bandwidth-bound regime and test weak-scaling behavior.
