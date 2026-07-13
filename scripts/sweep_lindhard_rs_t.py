#!/usr/bin/env python3

# ==========================================================================
# MOSAIQ-LINDHARD ENGINE
# Copyright (c) 2026 Bona Sapiens, Inc. All rights reserved.
# * This software is licensed under the MosaiQ-Lindhard Source Code License Agreement.
# Free for non-commercial personal and academic research use only.
# Commercial, governmental, and public institutional use requires a
# separate paid license. See LICENSE file for details.
# * Contact: kim.ingee@bonasapiens.com
# ==========================================================================

"""(r_s, T) sweep validation of reverse-Dedekind singularity excision (Fig. 7)."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from plot_common import OUTPUT_DIR, apply_pdf_rcparams, output_path, save_figure

DATA_FILE = OUTPUT_DIR / "output_lindhard_rs_t_sweep.dat"
OUTPUT_FIG = output_path("lindhard_rs_t_excision")
SWEEP_BIN_CANDIDATES = (
    Path(__file__).resolve().parent.parent / "simulator" / "build" / "tests" / "sweep_lindhard_rs_t",
    Path(__file__).resolve().parent.parent / "simulator" / "build" / "sweep_lindhard_rs_t",
)


def find_sweep_binary() -> Path | None:
    for candidate in SWEEP_BIN_CANDIDATES:
        if candidate.is_file():
            return candidate
    return None


def run_sweep(binary: Path) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        [str(binary)],
        cwd=str(Path(__file__).resolve().parent.parent),
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        sys.stderr.write(completed.stdout)
        sys.stderr.write(completed.stderr)
        raise SystemExit(f"sweep binary failed with code {completed.returncode}")
    print(completed.stdout.strip())


def load_sweep(path: Path) -> dict[str, np.ndarray]:
    if not path.is_file():
        raise SystemExit(f"Error: data file not found: {path}")

    table = np.loadtxt(path, comments="#")
    if table.ndim != 2 or table.shape[1] < 10:
        raise SystemExit(f"Error: expected ≥10 columns in {path}, got shape {table.shape}")

    return {
        "r_s": table[:, 0],
        "T_K": table[:, 1],
        "tau": table[:, 2],
        "re_engine": table[:, 3],
        "im_engine": table[:, 4],
        "re_t0": table[:, 5],
        "im_t0": table[:, 6],
        "rel_err_re": table[:, 7],
        "rel_err_im": table[:, 8],
        "finite": table[:, 9],
    }


def pivot_grid(
    r_s: np.ndarray,
    t_k: np.ndarray,
    values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rs_unique = np.unique(r_s)
    t_unique = np.unique(t_k)
    grid = np.full((len(t_unique), len(rs_unique)), np.nan, dtype=float)
    rs_index = {float(v): i for i, v in enumerate(rs_unique)}
    t_index = {float(v): i for i, v in enumerate(t_unique)}
    for rs_val, t_val, value in zip(r_s, t_k, values, strict=True):
        grid[t_index[float(t_val)], rs_index[float(rs_val)]] = value
    return rs_unique, t_unique, grid


def configure_matplotlib() -> None:
    apply_pdf_rcparams(
        {
            "font.family": "serif",
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
        }
    )


def render_heatmap(data: dict[str, np.ndarray], destination: Path) -> None:
    configure_matplotlib()

    rs_u, t_u, err_re = pivot_grid(data["r_s"], data["T_K"], data["rel_err_re"])
    _, _, finite = pivot_grid(data["r_s"], data["T_K"], data["finite"])

    # Floor for log color scale; zero error cells stay visible as the deep end.
    plot_err = np.where(np.isfinite(err_re), np.maximum(err_re, 1.0e-16), np.nan)
    nan_mask = ~np.isfinite(finite) | (finite < 0.5) | ~np.isfinite(err_re)

    fig, ax = plt.subplots(figsize=(8.2, 5.8), layout="constrained")
    rs_mesh, t_mesh = np.meshgrid(rs_u, t_u)

    finite_vals = plot_err[~nan_mask]
    if finite_vals.size == 0:
        raise SystemExit("No finite cells to plot — excision failed across the sweep.")

    norm = LogNorm(vmin=max(finite_vals.min(), 1.0e-16), vmax=max(finite_vals.max(), 1.0e-8))
    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad("#3d0000")

    mesh = ax.pcolormesh(
        rs_mesh,
        t_mesh,
        np.ma.masked_where(nan_mask, plot_err),
        shading="auto",
        cmap=cmap,
        norm=norm,
    )
    contours = ax.contour(
        rs_u,
        t_u,
        np.ma.masked_where(nan_mask, plot_err),
        levels=np.logspace(-12, -1, 8),
        colors="white",
        linewidths=0.45,
        alpha=0.55,
    )
    ax.clabel(contours, fmt=lambda v: f"{v:.0e}", fontsize=7, colors="white")

    ax.set_yscale("log")
    ax.set_xlabel(r"$r_s$")
    ax.set_ylabel(r"$T\,[\mathrm{K}]$")
    ax.set_title(
        r"Reverse-Dedekind singularity excision: "
        r"$|\Re\chi^L(T)-\Re\chi^L_{T=0}|/|\Re\chi^L_{T=0}|$"
        "\n"
        r"$(\bar{q},\bar{\omega})=(1.0,0.5),\ \gamma=1$"
    )
    ax.set_xlim(rs_u.min(), rs_u.max())
    ax.set_ylim(t_u.min(), t_u.max())

    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(r"relative error vs causality-pinned $T=0$")

    n_nan = int(np.count_nonzero(nan_mask))
    ax.text(
        0.02,
        0.02,
        f"NaN/Inf cells: {n_nan}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        color="white",
        bbox={"facecolor": "black", "alpha": 0.35, "pad": 3, "edgecolor": "none"},
    )

    save_figure(fig, destination)
    plt.close(fig)
    print(f"Wrote {destination}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--skip-run",
        action="store_true",
        help="Plot existing output/output_lindhard_rs_t_sweep.dat without re-running the C++ sweep",
    )
    args = parser.parse_args()

    if not args.skip_run:
        binary = find_sweep_binary()
        if binary is None:
            raise SystemExit(
                "sweep_lindhard_rs_t binary not found. Build simulator/tests first "
                "(cmake --build simulator/build --target sweep_lindhard_rs_t)."
            )
        run_sweep(binary)

    data = load_sweep(DATA_FILE)
    n_nan = int(np.count_nonzero(data["finite"] < 0.5))
    print(
        f"Loaded {DATA_FILE.name}: {data['r_s'].size} samples, "
        f"NaN/Inf={n_nan}, "
        f"max rel_err_re={np.nanmax(data['rel_err_re']):.3e}"
    )
    render_heatmap(data, OUTPUT_FIG)


if __name__ == "__main__":
    main()
