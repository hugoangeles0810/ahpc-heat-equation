#!/usr/bin/env python3
"""Compute and plot performance metrics for the MPI heat solver.

Reads the benchmark CSV produced by scripts/collect_results.sh and produces
three figures (objective c of the assignment):

  * execution time      vs. number of processes
  * speedup   S(p)=T(1)/T(p)   vs. p   (with the ideal linear reference)
  * efficiency E(p)=S(p)/p     vs. p

Usage:
    python analysis/plot_performance.py                 # defaults below
    python analysis/plot_performance.py --csv results/benchmark.csv \
                                        --outdir results

Expected CSV columns:
    nprocs,idim,kdim,iterations,wall_time_s,comm_time_s,criterion_time_s
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def load(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    # average repeated runs at the same process count, then sort by p
    df = df.groupby("nprocs", as_index=False).mean(numeric_only=True)
    df = df.sort_values("nprocs").reset_index(drop=True)

    baseline = df.loc[df["nprocs"].idxmin()]
    t1 = baseline["wall_time_s"]
    df["speedup"] = t1 / df["wall_time_s"]
    df["efficiency"] = df["speedup"] / df["nprocs"]
    return df


def plot_time(df: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(df["nprocs"], df["wall_time_s"], "o-", label="wall clock")
    if "comm_time_s" in df:
        ax.plot(df["nprocs"], df["comm_time_s"], "s--", label="communication")
    ax.set(xlabel="processes (p)", ylabel="time [s]", title="Execution time")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.grid(True, which="both", ls=":")
    ax.legend()
    return _save(fig, outdir / "time.png")


def plot_speedup(df: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(df["nprocs"], df["speedup"], "o-", label="measured")
    ax.plot(df["nprocs"], df["nprocs"], "k--", alpha=0.5, label="ideal (linear)")
    ax.set(xlabel="processes (p)", ylabel="speedup S(p)", title="Speedup")
    ax.grid(True, ls=":")
    ax.legend()
    return _save(fig, outdir / "speedup.png")


def plot_efficiency(df: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(df["nprocs"], df["efficiency"], "o-", label="measured")
    ax.axhline(1.0, color="k", ls="--", alpha=0.5, label="ideal (1.0)")
    ax.set(xlabel="processes (p)", ylabel="efficiency E(p)", title="Efficiency",
           ylim=(0, 1.1))
    ax.grid(True, ls=":")
    ax.legend()
    return _save(fig, outdir / "efficiency.png")


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

    print(df[["nprocs", "wall_time_s", "speedup", "efficiency"]].to_string(index=False))
    for fn in (plot_time, plot_speedup, plot_efficiency):
        print("wrote", fn(df, args.outdir))


if __name__ == "__main__":
    main()
