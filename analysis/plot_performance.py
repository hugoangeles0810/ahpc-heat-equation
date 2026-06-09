#!/usr/bin/env python3
"""Compute and plot performance metrics for the MPI heat solver.

Reads the benchmark CSV produced by scripts/collect_results.sh and produces
three figures (objective c of the assignment):

  * execution time      vs. number of processes
  * speedup   S(p)=T(1)/T(p)   vs. p   (with the ideal linear reference)
  * efficiency E(p)=S(p)/p     vs. p

The benchmark sweeps two axes (grid size N and process count p), so every
metric is computed *per grid size*: speedup/efficiency for a given n use that
n's own single-process run as the baseline, and each figure draws one curve
per grid size.

Usage:
    python analysis/plot_performance.py                 # defaults below
    python analysis/plot_performance.py --csv results/benchmark.csv \
                                        --outdir results

Expected CSV columns (header written once by scripts/benchmark.slurm):
    grid_size,nprocs,idim,kdim,iterations,wall_time_s,comm_time_s,criterion_time_s
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def load(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)

    # average repeated runs at the same (grid_size, nprocs), then sort
    df = df.groupby(["grid_size", "nprocs"], as_index=False).mean(numeric_only=True)
    df = df.sort_values(["grid_size", "nprocs"]).reset_index(drop=True)

    # speedup/efficiency are relative to the smallest p *within each grid size*.
    # The frame is sorted by nprocs, so "first" per grid is the baseline run.
    t1 = df.groupby("grid_size")["wall_time_s"].transform("first")
    df["speedup"] = t1 / df["wall_time_s"]
    df["efficiency"] = df["speedup"] / df["nprocs"]
    return df


def _grids(df: pd.DataFrame):
    """Yield (label, group) per grid size, sorted by n."""
    for n, g in df.groupby("grid_size"):
        yield f"n={n}", g


def plot_time(df: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(6, 4))
    for label, g in _grids(df):
        ax.plot(g["nprocs"], g["wall_time_s"], "o-", label=label)
    ax.set(xlabel="processes (p)", ylabel="time [s]", title="Execution time")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.grid(True, which="both", ls=":")
    ax.legend(title="grid size")
    return _save(fig, outdir / "time.png")


def plot_speedup(df: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(6, 4))
    for label, g in _grids(df):
        ax.plot(g["nprocs"], g["speedup"], "o-", label=label)
    procs = sorted(df["nprocs"].unique())
    ax.plot(procs, procs, "k--", alpha=0.5, label="ideal (linear)")
    ax.set(xlabel="processes (p)", ylabel="speedup S(p)", title="Speedup")
    ax.grid(True, ls=":")
    ax.legend(title="grid size")
    return _save(fig, outdir / "speedup.png")


def plot_efficiency(df: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(6, 4))
    for label, g in _grids(df):
        ax.plot(g["nprocs"], g["efficiency"], "o-", label=label)
    ax.axhline(1.0, color="k", ls="--", alpha=0.5, label="ideal (1.0)")
    ax.set(xlabel="processes (p)", ylabel="efficiency E(p)", title="Efficiency",
           ylim=(0, 1.1))
    ax.grid(True, ls=":")
    ax.legend(title="grid size")
    return _save(fig, outdir / "efficiency.png")


def plot_iterations(df: pd.DataFrame, outdir: Path) -> Path | None:
    """Iterations-to-convergence vs. grid size, validating the T=Θ(n²) bound.

    Convergence is a property of the numerics, not of p, so iterations are
    averaged over process counts for each grid size. A reference c·n² curve
    (anchored on the smallest grid) is overlaid; on log-log axes the data
    should track its slope of 2. Needs ≥2 distinct grid sizes to be meaningful.
    """
    by_n = df.groupby("grid_size", as_index=False)["iterations"].mean()
    by_n = by_n.sort_values(["grid_size"])
    if len(by_n) < 2:
        return None

    n = by_n["grid_size"]
    iters = by_n["iterations"]
    c = (iters / n**2).iloc[0]   # anchor c·n² on the smallest grid

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(n, iters, "o-", label="measured")
    ax.plot(n, c * n**2, "k--", alpha=0.5, label=r"$\Theta(n^2)$ reference")
    ax.set(xlabel="grid size (n)", ylabel="iterations to converge",
           title="Convergence cost vs. problem size")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.grid(True, which="both", ls=":")
    ax.legend()
    return _save(fig, outdir / "iterations.png")


def _save(fig, path: Path) -> Path:
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", type=Path, default=Path("results/benchmark.csv"))
    parser.add_argument("--outdir", type=Path, default=Path("results"))
    args = parser.parse_args()

    if not args.csv.exists():
        raise SystemExit(f"benchmark CSV not found: {args.csv}\n"
                         f"run scripts/benchmark.slurm first (or pass --csv).")

    args.outdir.mkdir(parents=True, exist_ok=True)
    df = load(args.csv)

    print(df[["grid_size", "nprocs", "wall_time_s", "speedup", "efficiency"]].to_string(index=False))
    for fn in (plot_time, plot_speedup, plot_efficiency, plot_iterations):
        path = fn(df, args.outdir)
        if path is None:
            print("skipped", fn.__name__, "(need ≥2 grid sizes)")
        else:
            print("wrote", path)


if __name__ == "__main__":
    main()
