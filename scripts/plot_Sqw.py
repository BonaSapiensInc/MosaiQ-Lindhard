#!/usr/bin/env python3
"""Render S_ee, S_ii, S_ei(q, omega) from SOPHUS multi-component RPA output."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import SymLogNorm

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_FILE = PROJECT_ROOT / "output" / "output_structure_factor.dat"
OUTPUT_DIR = PROJECT_ROOT / "output"

PLOT_SPECS: tuple[tuple[str, int, str], ...] = (
    ("S_ee", 3, r"$S_{ee}(q, \omega)$ — electron-electron"),
    ("S_ii", 5, r"$S_{ii}(q, \omega)$ — ion-ion (ion-acoustic)"),
    ("S_ei", 7, r"$S_{ei}(q, \omega)$ — electron-ion cross"),
)


def load_gridded(data: np.ndarray, value_col: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    q_1d = data[:, 0]
    w_1d = data[:, 1]
    values = data[:, value_col]

    q_unique = np.unique(q_1d)
    w_unique = np.unique(w_1d)
    q_grid, w_grid = np.meshgrid(q_unique, w_unique)
    value_grid = values.reshape(len(q_unique), len(w_unique)).T
    return q_grid, w_grid, value_grid


def choose_norm(grid: np.ndarray) -> SymLogNorm | None:
    finite = grid[np.isfinite(grid)]
    if finite.size == 0:
        return None

    abs_max = float(np.max(np.abs(finite)))
    if abs_max == 0.0:
        return None

    positive = finite[finite > 0.0]
    min_pos = float(np.min(positive)) if positive.size else abs_max * 1.0e-6
    linthresh = max(min_pos, abs_max * 1.0e-4, 1.0e-12)
    return SymLogNorm(linthresh=linthresh, linscale=1.0, vmin=-abs_max, vmax=abs_max)


def render_component(
    q_grid: np.ndarray,
    w_grid: np.ndarray,
    grid: np.ndarray,
    label: str,
    title: str,
    output_path: Path,
) -> None:
    norm = choose_norm(grid)
    plot_grid = np.ma.masked_invalid(grid)

    plt.rcParams.update({"font.size": 12, "font.family": "serif"})
    fig = plt.figure(figsize=(12, 10))

    ax1 = fig.add_subplot(211)
    contour_kwargs: dict = {"levels": 150, "cmap": cm.magma}
    if norm is not None:
        contour_kwargs["norm"] = norm
    contour = ax1.contourf(q_grid, w_grid, plot_grid, **contour_kwargs)
    cbar1 = fig.colorbar(contour, ax=ax1)
    cbar1.set_label(label)
    ax1.set_xlabel(r"Wave vector $q \ [k_F]$")
    ax1.set_ylabel(r"Frequency $\omega \ [E_F/\hbar]$")
    ax1.set_title(f"{title} — Contour Map", fontweight="bold")

    ax2 = fig.add_subplot(212, projection="3d")
    z_plot = np.nan_to_num(plot_grid, nan=0.0)
    if norm is not None:
        z_plot = np.clip(z_plot, norm.vmin, norm.vmax)
    surface = ax2.plot_surface(
        q_grid,
        w_grid,
        z_plot,
        cmap=cm.magma,
        norm=norm,
        edgecolor="none",
        antialiased=True,
        alpha=0.92,
    )
    cbar2 = fig.colorbar(surface, ax=ax2, pad=0.1)
    cbar2.set_label(label)
    ax2.set_xlabel(r"$q \ [k_F]$", labelpad=10)
    ax2.set_ylabel(r"$\omega \ [E_F/\hbar]$", labelpad=10)
    ax2.set_zlabel(label, labelpad=10)
    ax2.set_title(f"{title} — 3D Surface", fontweight="bold")
    ax2.view_init(elev=35, azim=-120)

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output_path}")


def main() -> None:
    if not DATA_FILE.is_file():
        raise SystemExit(f"Error: data file not found: {DATA_FILE}")

    data = np.loadtxt(DATA_FILE, comments="#")
    if data.ndim != 2 or data.shape[1] < 8:
        raise SystemExit(
            f"Error: expected at least 8 columns in {DATA_FILE}, got shape {data.shape}"
        )

    for stem, col, title in PLOT_SPECS:
        q_grid, w_grid, value_grid = load_gridded(data, col)
        render_component(
            q_grid,
            w_grid,
            value_grid,
            label=title.split("—")[0].strip(),
            title=title,
            output_path=OUTPUT_DIR / f"{stem}_plot.png",
        )


if __name__ == "__main__":
    main()
