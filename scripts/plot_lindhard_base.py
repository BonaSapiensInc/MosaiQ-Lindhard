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

"""Render bare Lindhard response manifolds (Figure 3) as six contour panels."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from plot_common import OUTPUT_DIR, PDF_RCPARAMS, output_path
from plot_Sqw import (
    ION_CHANNEL_OMEGA_MAX,
    ION_CHANNEL_OMEGA_MIN,
    configure_seaborn,
    load_gridded,
    prepare_plot_grid,
    render_contour,
)

DATA_FILE = OUTPUT_DIR / "output_lindhard_base.dat"

ELECTRON_Q_MAX = 4.0
ION_Q_MAX = 50.0
ELECTRON_OMEGA_MAX = 3.0

# (output stem, value column, q cap, colorbar label, title)
PANEL_SPECS: tuple[tuple[str, int, float, str, str], ...] = (
    (
        "chi_e_re",
        2,
        ELECTRON_Q_MAX,
        r"$|\Re \tilde{\chi}_e^L|$",
        r"$\Re \tilde{\chi}_e^L(\bar{q}, \bar{\omega})$ — electron",
    ),
    (
        "chi_e_im",
        3,
        ELECTRON_Q_MAX,
        r"$|\Im \tilde{\chi}_e^L|$",
        r"$\Im \tilde{\chi}_e^L(\bar{q}, \bar{\omega})$ — electron",
    ),
    (
        "chi_i_re",
        4,
        ION_Q_MAX,
        r"$|\Re \tilde{\chi}_i^L|$",
        r"$\Re \tilde{\chi}_i^L(\bar{q}, \bar{\omega})$ — ion",
    ),
    (
        "chi_i_im",
        5,
        ION_Q_MAX,
        r"$|\Im \tilde{\chi}_i^L|$",
        r"$\Im \tilde{\chi}_i^L(\bar{q}, \bar{\omega})$ — ion",
    ),
    (
        "chi_ei_re",
        6,
        ELECTRON_Q_MAX,
        r"$|\Re \tilde{\chi}_{ei}^L|$",
        r"$\Re \tilde{\chi}_{ei}^L(\bar{q}, \bar{\omega})$ — cross",
    ),
    (
        "chi_ei_im",
        7,
        ELECTRON_Q_MAX,
        r"$|\Im \tilde{\chi}_{ei}^L|$",
        r"$\Im \tilde{\chi}_{ei}^L(\bar{q}, \bar{\omega})$ — cross",
    ),
)


def load_data(path: Path) -> np.ndarray:
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing {path}. Run the MosaiQ simulator (e.g. ./simulator/build/mosaiq 1.0 10000)."
        )
    return np.loadtxt(path, comments="#")


def channel_grid(
    data: np.ndarray,
    value_col: int,
    q_max: float,
    omega_max: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = data[data[:, 0] <= q_max + 1.0e-9]
    if omega_max is not None:
        rows = rows[rows[:, 1] <= omega_max + 1.0e-9]
    q_grid, w_grid, value_grid = load_gridded(rows, value_col)
    return q_grid, w_grid, np.abs(value_grid)


def main() -> None:
    configure_seaborn()
    data = load_data(DATA_FILE)

    for stem, value_col, q_max, cbar_label, title in PANEL_SPECS:
        ion_panel = stem.startswith("chi_i")
        omega_cap = ION_CHANNEL_OMEGA_MAX if ion_panel else ELECTRON_OMEGA_MAX
        q_grid, w_grid, magnitude = channel_grid(
            data, value_col, q_max, omega_max=omega_cap if ion_panel else None
        )
        plot_grid, color_norm = prepare_plot_grid(magnitude)
        if color_norm is None:
            raise RuntimeError(f"No finite positive amplitudes for {stem}")

        render_contour(
            q_grid,
            w_grid,
            plot_grid,
            color_norm,
            cbar_label,
            title,
            output_path(stem),
            q_xmax=q_max,
            omega_ymin=ION_CHANNEL_OMEGA_MIN if ion_panel else None,
            omega_ymax=omega_cap,
            contour_profile="cross" if stem.startswith("chi_ei") else "default",
        )


if __name__ == "__main__":
    main()
