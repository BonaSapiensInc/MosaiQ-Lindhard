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

import numpy as np

from plot_common import OUTPUT_DIR, output_path
from plot_Sqw import (
    ELECTRON_DISPLAY_BOUNDS,
    ION_DISPLAY_BOUNDS,
    configure_seaborn,
    interpolate_to_display_mesh,
    load_gridded,
    prepare_ion_channel_mesh,
    prepare_ion_plot_grid,
    prepare_plot_grid,
    render_contour,
)

DATA_FILE = OUTPUT_DIR / "output_lindhard_base.dat"

ELECTRON_Q_MAX = 4.0
ION_Q_MAX = 50.0

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


def main() -> None:
    configure_seaborn()
    data = load_data(DATA_FILE)

    for stem, value_col, q_max, cbar_label, title in PANEL_SPECS:
        ion_panel = stem.startswith("chi_i")
        rows = data[data[:, 0] <= q_max + 1.0e-9]
        q_grid, w_grid, magnitude = load_gridded(rows, value_col)
        magnitude = np.abs(magnitude)

        if ion_panel:
            q_display, w_display, value_display = prepare_ion_channel_mesh(
                q_grid, w_grid, magnitude
            )
            plot_grid, color_norm = prepare_ion_plot_grid(
                q_display, w_display, value_display
            )
            q_bounds = ION_DISPLAY_BOUNDS[:2]
            w_bounds = ION_DISPLAY_BOUNDS[2:]
        else:
            q_display, w_display, value_display = interpolate_to_display_mesh(
                q_grid, w_grid, magnitude, *ELECTRON_DISPLAY_BOUNDS
            )
            plot_grid, color_norm = prepare_plot_grid(value_display)
            q_bounds = ELECTRON_DISPLAY_BOUNDS[:2]
            w_bounds = ELECTRON_DISPLAY_BOUNDS[2:]

        if color_norm is None:
            raise RuntimeError(f"No finite positive amplitudes for {stem}")

        render_contour(
            q_display,
            w_display,
            plot_grid,
            color_norm,
            cbar_label,
            title,
            output_path(stem),
            q_bounds=q_bounds,
            w_bounds=w_bounds,
            contour_profile="cross" if stem.startswith("chi_ei") else "default",
        )


if __name__ == "__main__":
    main()
