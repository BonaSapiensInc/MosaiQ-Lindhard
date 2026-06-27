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

"""Publication-quality plasmon dispersion and Landau damping figure from MosaiQ-Lindhard output (PRE style)."""

from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib import cm
from matplotlib.colors import LogNorm, Normalize
from matplotlib.figure import Figure
import numpy as np

from plot_common import OUTPUT_DIR, apply_pdf_rcparams, output_path, save_figure
from plot_Sqw import build_contour_grid, log_contour_levels, prepare_plot_grid

DISPERSION_FILE = OUTPUT_DIR / "output_plasmon_dispersion.dat"
STRUCTURE_FACTOR_FILE = OUTPUT_DIR / "output_structure_factor.dat"
OUTPUT_FIG = output_path("plasmon_dispersion")

S_EE_COLUMN = 3

PANEL_A_OMEGA_YMAX = 10.0

COLOR_OMEGA_P = "#111111"
COLOR_BOHM_GROSS = "#ffffff"
COLOR_BOHM_GROSS_LEGEND = "#555555"
COLOR_LANDAU = "#8b1e1e"
COLOR_SPE_CHANNEL = "#2f5f8f"
COLOR_SPE_LINE = "#174a73"
SPE_FILL_ALPHA = 0.30
LOG_FLOOR = 1.0e-12


def bohm_gross_path_effects() -> list[pe.AbstractPathEffect]:
    """Dark halo so the white dashed curve reads on rainbow backgrounds and in the legend."""
    return [
        pe.Stroke(linewidth=3.6, foreground="#2b2b2b"),
        pe.Normal(),
    ]


def bohm_gross_legend_handle() -> plt.Line2D:
    line = plt.Line2D(
        [0],
        [0],
        color=COLOR_BOHM_GROSS_LEGEND,
        linewidth=2.0,
        linestyle="--",
        dashes=(6, 3),
    )
    return line


def load_dispersion(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if not path.is_file():
        raise SystemExit(f"Error: data file not found: {path}")

    table = np.loadtxt(path, comments="#")
    if table.ndim != 2 or table.shape[1] < 4:
        raise SystemExit(f"Error: expected 4 columns in {path}, got shape {table.shape}")

    q = table[:, 0]
    omega_p = table[:, 1]
    landau = table[:, 2]
    bohm_gross = table[:, 3]
    return q, omega_p, landau, bohm_gross


def parse_electron_reduced_params(path: Path) -> tuple[float, float]:
    pattern = re.compile(r"gamma_e\s*=\s*([^\s]+)\s+tau_e\s*=\s*([^\s]+)")
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            match = pattern.search(line)
            if match:
                return float(match.group(1)), float(match.group(2))
    raise SystemExit(f"Error: could not parse gamma_e/tau_e from {path}")


def bare_imaginary_lindhard_along_roots(
    q: np.ndarray,
    omega_p: np.ndarray,
    tau_e: float,
    gamma_e: float,
) -> np.ndarray:
    """Magnitude of bare Im[chi_e^L] on the extracted plasmon branch (dissipation weight)."""
    bare = np.zeros_like(q, dtype=float)
    mask = np.isfinite(omega_p) & (q > 0.0)
    if not np.any(mask):
        return bare

    q_vals = q[mask]
    omega_vals = omega_p[mask]
    nu_p = ((omega_vals / q_vals) + q_vals) / 2.0
    nu_n = ((omega_vals / q_vals) - q_vals) / 2.0

    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        log_ratio = np.log((1.0 + np.exp(-(nu_n * nu_n - gamma_e) / tau_e)) /
                           (1.0 + np.exp(-(nu_p * nu_p - gamma_e) / tau_e)))
        bare[mask] = -(np.pi / q_vals) * tau_e * log_ratio

    return np.abs(bare)


def load_s_ee_background(
    path: Path,
    q: np.ndarray,
    omega_p: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ma.MaskedArray, LogNorm | Normalize, float, float] | None:
    if not path.is_file():
        return None

    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] <= S_EE_COLUMN:
        return None

    active = np.isfinite(omega_p) & (q > 0.0)
    if np.any(active):
        q_bg_max = min(float(np.max(q[active]) * 1.20), float(np.max(data[:, 0])))
        w_bg_max = PANEL_A_OMEGA_YMAX
    else:
        q_bg_max = min(4.0, float(np.max(data[:, 0])))
        w_bg_max = PANEL_A_OMEGA_YMAX

    q_grid, w_grid, value_grid = build_contour_grid(
        data, S_EE_COLUMN, 0.1, q_bg_max, w_bg_max, n_omega=250
    )
    plot_grid, norm = prepare_plot_grid(value_grid)
    color_norm: LogNorm | Normalize = norm if norm is not None else Normalize()
    return q_grid, w_grid, plot_grid, color_norm, q_bg_max, w_bg_max


def finite_segments(q: np.ndarray, y: np.ndarray) -> list[tuple[np.ndarray, np.ndarray]]:
    """Split (q, y) into contiguous segments where y is finite."""
    mask = np.isfinite(y)
    segments: list[tuple[np.ndarray, np.ndarray]] = []
    start: int | None = None

    for index, is_finite in enumerate(mask):
        if is_finite and start is None:
            start = index
        elif not is_finite and start is not None:
            segments.append((q[start:index], y[start:index]))
            start = None

    if start is not None:
        segments.append((q[start:], y[start:]))

    return segments


def configure_matplotlib() -> None:
    apply_pdf_rcparams(
        {
            "font.family": "serif",
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "legend.fontsize": 10,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.linewidth": 0.8,
            "grid.alpha": 0.25,
            "grid.linewidth": 0.5,
        }
    )


def render_s_ee_background(
    ax: plt.Axes,
    q_grid: np.ndarray,
    w_grid: np.ndarray,
    plot_grid: np.ma.MaskedArray,
    color_norm: LogNorm | Normalize,
) -> None:
    cmap = cm.rainbow
    if isinstance(color_norm, LogNorm):
        fill_levels, line_levels = log_contour_levels(color_norm, n_filled=200, n_lines=48)
        ax.contourf(
            q_grid,
            w_grid,
            plot_grid,
            levels=fill_levels,
            cmap=cmap,
            norm=color_norm,
            zorder=0,
        )
        ax.contour(
            q_grid,
            w_grid,
            plot_grid,
            levels=line_levels,
            norm=color_norm,
            colors="white",
            linewidths=0.35,
            alpha=0.65,
            zorder=1,
        )
    else:
        ax.contourf(
            q_grid,
            w_grid,
            plot_grid,
            levels=150,
            cmap=cmap,
            norm=color_norm,
            zorder=0,
        )


def plot_dispersion(
    q: np.ndarray,
    omega_p: np.ndarray,
    landau: np.ndarray,
    bohm_gross: np.ndarray,
    background: tuple[np.ndarray, np.ndarray, np.ma.MaskedArray, LogNorm | Normalize, float, float] | None,
    bare_spe: np.ndarray,
) -> Figure:
    configure_matplotlib()

    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=(7.0, 5.8),
        gridspec_kw={"height_ratios": [1.45, 1.0], "hspace": 0.08},
    )

    q_view_max = float(np.max(q[np.isfinite(omega_p)])) * 1.12 if np.any(np.isfinite(omega_p)) else 4.0

    if background is not None:
        q_grid, w_grid, plot_grid, color_norm, q_bg_max, w_bg_max = background
        render_s_ee_background(ax_top, q_grid, w_grid, plot_grid, color_norm)
        ax_top.set_facecolor("0.15")
        q_view_max = max(q_view_max, q_bg_max)

    for q_seg, y_seg in finite_segments(q, bohm_gross):
        (line,) = ax_top.plot(
            q_seg,
            y_seg,
            color=COLOR_BOHM_GROSS,
            linewidth=2.0,
            linestyle="--",
            dashes=(6, 3),
            solid_capstyle="round",
            zorder=4,
        )
        line.set_path_effects(bohm_gross_path_effects())

    for q_seg, y_seg in finite_segments(q, omega_p):
        ax_top.plot(
            q_seg,
            y_seg,
            color=COLOR_OMEGA_P,
            linewidth=2.2,
            solid_capstyle="round",
            zorder=5,
        )

    ax_top.set_ylabel(r"$\bar{\omega}_p = \hbar\omega_p/\epsilon_\mathrm{F}$")
    ax_top.set_ylim(0.0, PANEL_A_OMEGA_YMAX)
    ax_top.set_xlim(0.0, q_view_max)
    ax_top.grid(False)

    legend_handles = [
        plt.Line2D([0], [0], color=COLOR_OMEGA_P, linewidth=2.2),
        bohm_gross_legend_handle(),
    ]
    legend_labels = [r"Extracted $\bar{\omega}_p$", "Bohm--Gross"]
    if background is not None:
        legend_handles.insert(
            0,
            plt.Rectangle((0, 0), 1, 1, facecolor=cm.rainbow(0.65), edgecolor="0.4"),
        )
        legend_labels.insert(0, r"$S_{ee}(\bar{q}, \bar{\omega})$")

    ax_top.legend(
        legend_handles,
        legend_labels,
        loc="upper left",
        frameon=True,
        framealpha=0.92,
        edgecolor="0.75",
        facecolor="white",
    )

    spe_active = np.isfinite(omega_p) & np.isfinite(bare_spe) & (bare_spe > 0.0)
    ax_spe = ax_bottom.twinx()
    if np.any(spe_active):
        spe_curve = np.where(spe_active, bare_spe, np.nan)
        ax_spe.fill_between(
            q,
            0.0,
            spe_curve,
            where=spe_active,
            color=COLOR_SPE_CHANNEL,
            alpha=SPE_FILL_ALPHA,
            interpolate=True,
            zorder=0,
        )
        for q_seg, y_seg in finite_segments(q, spe_curve):
            ax_spe.plot(
                q_seg,
                y_seg,
                color=COLOR_SPE_LINE,
                linewidth=2.0,
                linestyle="-",
                solid_capstyle="round",
                zorder=1,
            )

        spe_max = float(np.nanmax(bare_spe[spe_active]))
        ax_spe.set_ylim(0.0, max(spe_max * 1.12, 0.5))
        ax_spe.set_ylabel(
            r"$|\Im\chi_e^L(\bar{q}, \bar{\omega}_p)|$",
            color=COLOR_SPE_LINE,
        )
        ax_spe.tick_params(axis="y", colors=COLOR_SPE_LINE)
        ax_spe.spines["right"].set_color(COLOR_SPE_LINE)

    for q_seg, y_seg in finite_segments(q, landau):
        positive = y_seg > 0.0
        if np.count_nonzero(positive) >= 2:
            ax_bottom.plot(
                q_seg[positive],
                y_seg[positive],
                color=COLOR_LANDAU,
                linewidth=2.0,
                solid_capstyle="round",
                zorder=3,
            )

    ax_bottom.set_yscale("log")
    ax_bottom.set_ylabel(r"$\Im[\varepsilon(\bar{q}, \bar{\omega}_p)]$", color=COLOR_LANDAU)
    ax_bottom.tick_params(axis="y", colors=COLOR_LANDAU)
    ax_bottom.spines["left"].set_color(COLOR_LANDAU)
    ax_bottom.set_xlabel(r"Reduced wavevector $\bar{q} = q/k_\mathrm{F}$")
    ax_bottom.set_xlim(0.0, q_view_max)
    ax_bottom.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.35)

    bottom_handles = [
        plt.Rectangle(
            (0, 0),
            1,
            1,
            facecolor=COLOR_SPE_CHANNEL,
            alpha=SPE_FILL_ALPHA,
            edgecolor=COLOR_SPE_LINE,
        ),
        plt.Line2D([0], [0], color=COLOR_SPE_LINE, linewidth=2.0),
        plt.Line2D([0], [0], color=COLOR_LANDAU, linewidth=2.0),
    ]
    bottom_labels = [
        r"Bare SPE channel $|\Im\chi_e^L(\bar{q}, \bar{\omega}_p)|$",
        r"SPE weight on $\bar{\omega}_p$",
        r"$\Im[\varepsilon(\bar{q}, \bar{\omega}_p)]$",
    ]
    ax_bottom.legend(
        bottom_handles,
        bottom_labels,
        loc="upper right",
        frameon=True,
        framealpha=0.94,
        edgecolor="0.75",
        facecolor="white",
        fontsize=9,
    )

    ax_top.text(
        0.98,
        0.06,
        "(a)",
        transform=ax_top.transAxes,
        color="white",
        fontsize=12,
        fontweight="bold",
        va="bottom",
        ha="right",
        zorder=10,
    )
    ax_bottom.text(
        0.98,
        0.08,
        "(b)",
        transform=ax_bottom.transAxes,
        color="black",
        fontsize=12,
        fontweight="bold",
        va="bottom",
        ha="right",
        zorder=10,
    )

    fig.align_ylabels([ax_top, ax_bottom])
    fig.subplots_adjust(hspace=0.08)
    return fig


def main() -> None:
    q, omega_p, landau, bohm_gross = load_dispersion(DISPERSION_FILE)
    background = load_s_ee_background(STRUCTURE_FACTOR_FILE, q, omega_p)

    params_path = DISPERSION_FILE if DISPERSION_FILE.is_file() else STRUCTURE_FACTOR_FILE
    gamma_e, tau_e = parse_electron_reduced_params(params_path)
    bare_spe = bare_imaginary_lindhard_along_roots(q, omega_p, tau_e, gamma_e)

    fig = plot_dispersion(q, omega_p, landau, bohm_gross, background, bare_spe)

    saved_path = save_figure(fig, OUTPUT_FIG.stem)
    plt.close(fig)

    n_roots = int(np.sum(np.isfinite(omega_p)))
    print(f"Loaded {q.size} q points ({n_roots} finite plasmon roots).")
    if background is None:
        print(f"Warning: no background map from {STRUCTURE_FACTOR_FILE}; dispersion-only plot saved.")
    print(f"Saved {saved_path}")


if __name__ == "__main__":
    main()
