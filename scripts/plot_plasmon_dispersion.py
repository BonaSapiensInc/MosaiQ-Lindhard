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
from matplotlib import cm
from matplotlib.colors import LogNorm, Normalize
from matplotlib.figure import Figure
import numpy as np

from plot_common import OUTPUT_DIR, apply_pdf_rcparams, format_temperature_k, output_path, save_figure
from plot_Sqw import (
    CONTOUR_ISOLINE_ALPHA,
    CONTOUR_ISOLINE_LINEWIDTH,
    build_contour_grid,
    draw_contour_isolines,
    filter_contour_line_levels,
    frame_contour_axis,
    interpolate_to_display_mesh,
    log_contour_levels,
    prepare_plot_grid,
    rainbow_log_cmap,
)

# Manuscript Figure 5 parameters (must match between dispersion and structure-factor exports).
MANUSCRIPT_FIG5_RS = 2.0
MANUSCRIPT_FIG5_T_KELVIN = 10_000.0

OUTPUT_FIG = output_path("plasmon_dispersion")

S_EE_COLUMN = 3

PANEL_A_OMEGA_YMAX = 6.0
PANEL_A_Q_SNAP = 0.5
PANEL_A_OMEGA_SNAP = 0.5

RS_HEADER = re.compile(r"r_s\s*=\s*([^\s]+)")
T_HEADER = re.compile(r"T_K\s*=\s*([^\s]+)")

COLOR_OMEGA_P = "#111111"
COLOR_BOHM_GROSS = "#111111"
COLOR_LANDAU = "#8b1e1e"
COLOR_SPE_CHANNEL = "#2f5f8f"
COLOR_SPE_LINE = "#174a73"
SPE_FILL_ALPHA = 0.30
LOG_FLOOR = 1.0e-12


def bohm_gross_legend_handle() -> plt.Line2D:
    return plt.Line2D(
        [0],
        [0],
        color=COLOR_BOHM_GROSS,
        linewidth=2.0,
        linestyle="--",
    )


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


def parse_run_metadata(path: Path) -> tuple[float | None, float | None]:
    rs_value: float | None = None
    t_value: float | None = None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            rs_match = RS_HEADER.search(line)
            if rs_match:
                rs_value = float(rs_match.group(1))
            t_match = T_HEADER.search(line)
            if t_match:
                t_value = float(t_match.group(1))
    return rs_value, t_value


def rs_filename_tag(rs: float) -> str:
    return str(int(rs)) if abs(rs - round(rs)) < 1.0e-9 else str(rs).replace(".", "p")


def resolve_fig5_data_paths(
    rs: float = MANUSCRIPT_FIG5_RS,
    t_kelvin: float = MANUSCRIPT_FIG5_T_KELVIN,
) -> tuple[Path, Path]:
    """Require dispersion and structure-factor exports from the same simulator run."""
    candidates: list[tuple[Path, Path]] = [
        (
            OUTPUT_DIR / f"output_plasmon_dispersion_rs{rs_filename_tag(rs)}.dat",
            OUTPUT_DIR / f"output_structure_factor_rs{rs_filename_tag(rs)}.dat",
        ),
        (
            OUTPUT_DIR / "output_plasmon_dispersion.dat",
            OUTPUT_DIR / "output_structure_factor.dat",
        ),
    ]

    for dispersion_path, structure_factor_path in candidates:
        if not dispersion_path.is_file() or not structure_factor_path.is_file():
            continue
        disp_rs, disp_t = parse_run_metadata(dispersion_path)
        sf_rs, sf_t = parse_run_metadata(structure_factor_path)
        if disp_rs is None or sf_rs is None:
            continue
        if abs(disp_rs - sf_rs) > 1.0e-6:
            continue
        if abs(disp_rs - rs) > 1.0e-6:
            continue
        if disp_t is not None and sf_t is not None and abs(disp_t - sf_t) > 1.0e-3:
            continue
        if disp_t is not None and abs(disp_t - t_kelvin) > 1.0e-3:
            continue
        return dispersion_path, structure_factor_path

    raise SystemExit(
        "Error: matched Figure 5 data not found.\n"
        f"  Expected r_s={rs}, T={t_kelvin} K in both:\n"
        "    output/output_plasmon_dispersion.dat\n"
        "    output/output_structure_factor.dat\n"
        "  Regenerate with:\n"
        f"    ./simulator/build/mosaiq_simulator {rs} {t_kelvin}\n"
        f"    cp output/output_plasmon_dispersion.dat output/output_plasmon_dispersion_rs{rs_filename_tag(rs)}.dat\n"
        f"    cp output/output_structure_factor.dat output/output_structure_factor_rs{rs_filename_tag(rs)}.dat"
    )


def snap_ceiling(value: float, step: float) -> float:
    return float(np.ceil(value / step) * step)


def compute_panel_a_bounds(
    q: np.ndarray,
    omega_p: np.ndarray,
) -> tuple[float, float, float, float]:
    """Tight display rectangle snapped to plasmon extent (fills axis edges)."""
    active = np.isfinite(omega_p) & (q > 0.0)
    if not np.any(active):
        return 0.0, 2.0, 0.0, PANEL_A_OMEGA_YMAX

    q_hi = snap_ceiling(float(np.max(q[active])), PANEL_A_Q_SNAP)
    omega_top = float(np.max(omega_p[active]))
    w_hi = max(PANEL_A_OMEGA_YMAX, snap_ceiling(omega_top, PANEL_A_OMEGA_SNAP))
    return 0.0, q_hi, 0.0, w_hi


def fill_display_boundary(grid: np.ndarray) -> np.ndarray:
    """Extend nearest valid samples to the display rectangle edges."""
    out = np.array(grid, dtype=float, copy=True)
    for i in range(out.shape[0]):
        row = out[i]
        finite = np.isfinite(row)
        if not np.any(finite):
            continue
        j_first = int(np.argmax(finite))
        j_last = len(row) - 1 - int(np.argmax(finite[::-1]))
        row[:j_first] = row[j_first]
        row[j_last + 1 :] = row[j_last]

    for j in range(out.shape[1]):
        col = out[:, j]
        finite = np.isfinite(col)
        if not np.any(finite):
            continue
        i_first = int(np.argmax(finite))
        i_last = len(col) - 1 - int(np.argmax(finite[::-1]))
        col[:i_first] = col[i_first]
        col[i_last + 1 :] = col[i_last]
    return out


def max_omega_in_structure_factor(path: Path) -> float:
    table = np.loadtxt(path, comments="#")
    return float(np.max(table[:, 1]))


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
    display_bounds: tuple[float, float, float, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, LogNorm | Normalize, float, float] | None:
    if not path.is_file():
        return None

    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] <= S_EE_COLUMN:
        return None

    q_lo, q_hi, w_lo, w_hi = display_bounds
    q_data_max = min(float(np.max(data[:, 0])), q_hi)
    if q_data_max <= 0.0:
        return None

    q_grid, w_grid, value_grid = build_contour_grid(
        data,
        S_EE_COLUMN,
        max(0.1, q_lo if q_lo > 0.0 else 0.1),
        q_data_max,
        w_hi,
        n_omega=400,
        w_min=max(w_lo, 0.01),
    )
    q_display, w_display, value_display = interpolate_to_display_mesh(
        q_grid,
        w_grid,
        np.abs(value_grid),
        q_lo,
        q_hi,
        w_lo,
        w_hi,
    )
    value_display = fill_display_boundary(value_display)
    plot_grid, norm = prepare_plot_grid(value_display)
    color_norm: LogNorm | Normalize = norm if norm is not None else Normalize()
    return q_display, w_display, plot_grid, color_norm, q_hi, w_hi


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
    plot_grid: np.ndarray,
    color_norm: LogNorm | Normalize,
) -> None:
    cmap = rainbow_log_cmap()
    if isinstance(color_norm, LogNorm):
        fill_levels, line_levels = log_contour_levels(color_norm, n_filled=200, n_lines=72)
        line_levels = filter_contour_line_levels(line_levels, color_norm)
        ax.contourf(
            q_grid,
            w_grid,
            plot_grid,
            levels=fill_levels,
            cmap=cmap,
            norm=color_norm,
            extend="min",
            zorder=0,
        )
        draw_contour_isolines(
            ax,
            q_grid,
            w_grid,
            plot_grid,
            color_norm,
            line_levels,
            line_width=CONTOUR_ISOLINE_LINEWIDTH,
            line_alpha=CONTOUR_ISOLINE_ALPHA,
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

    frame_contour_axis(ax)


def plot_dispersion(
    q: np.ndarray,
    omega_p: np.ndarray,
    landau: np.ndarray,
    bohm_gross: np.ndarray,
    background: tuple[np.ndarray, np.ndarray, np.ndarray, LogNorm | Normalize, float, float] | None,
    bare_spe: np.ndarray,
    display_bounds: tuple[float, float, float, float],
    temperature_k: float = MANUSCRIPT_FIG5_T_KELVIN,
) -> Figure:
    configure_matplotlib()

    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=(7.0, 5.8),
        gridspec_kw={"height_ratios": [1.45, 1.0], "hspace": 0.08},
    )

    q_lo, q_hi, w_lo, w_hi = display_bounds

    if background is not None:
        q_grid, w_grid, plot_grid, color_norm, _, _ = background
        render_s_ee_background(ax_top, q_grid, w_grid, plot_grid, color_norm)
        ax_top.set_xlim(q_lo, q_hi)
        ax_top.set_ylim(w_lo, w_hi)
        ax_top.margins(x=0, y=0)

    q_mask = (q >= q_lo) & (q <= q_hi)

    for q_seg, y_seg in finite_segments(q[q_mask], bohm_gross[q_mask]):
        ax_top.plot(
            q_seg,
            y_seg,
            color=COLOR_BOHM_GROSS,
            linewidth=2.0,
            linestyle="--",
            solid_capstyle="butt",
            dash_capstyle="butt",
            zorder=4,
        )

    for q_seg, y_seg in finite_segments(q[q_mask], omega_p[q_mask]):
        ax_top.plot(
            q_seg,
            y_seg,
            color=COLOR_OMEGA_P,
            linewidth=2.2,
            solid_capstyle="round",
            zorder=5,
        )

    ax_top.set_ylabel(r"$\bar{\omega}_p = \hbar\omega_p/\epsilon_\mathrm{F}$")
    if background is None:
        ax_top.set_ylim(w_lo, w_hi)
        ax_top.set_xlim(q_lo, q_hi)
        ax_top.margins(x=0, y=0)
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

    legend = ax_top.legend(
        legend_handles,
        legend_labels,
        loc="upper left",
        frameon=True,
        framealpha=0.92,
        edgecolor="0.75",
        facecolor="white",
    )

    fig.canvas.draw()
    legend_bbox = legend.get_window_extent(fig.canvas.get_renderer()).transformed(
        ax_top.transAxes.inverted()
    )
    ax_top.text(
        legend_bbox.x0,
        legend_bbox.y0 - 0.012,
        format_temperature_k(temperature_k),
        transform=ax_top.transAxes,
        color="white",
        fontsize=10,
        va="top",
        ha="left",
        zorder=10,
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
    ax_bottom.set_xlim(q_lo, q_hi)
    ax_bottom.margins(x=0)
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
        0.98,
        "(a)",
        transform=ax_top.transAxes,
        color="white",
        fontsize=12,
        fontweight="bold",
        va="top",
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
    dispersion_path, structure_factor_path = resolve_fig5_data_paths()
    omega_max = max_omega_in_structure_factor(structure_factor_path)
    if omega_max + 1.0e-6 < PANEL_A_OMEGA_YMAX:
        raise SystemExit(
            f"Error: {structure_factor_path} only reaches omega={omega_max:.3f}, "
            f"but panel (a) requires at least {PANEL_A_OMEGA_YMAX}.\n"
            "  Rebuild mosaiq_simulator and rerun the Figure 5 export."
        )

    q, omega_p, landau, bohm_gross = load_dispersion(dispersion_path)
    display_bounds = compute_panel_a_bounds(q, omega_p)
    _, _, _, w_hi = display_bounds
    if omega_max + 1.0e-6 < w_hi:
        raise SystemExit(
            f"Error: {structure_factor_path} only reaches omega={omega_max:.3f}, "
            f"but panel (a) requires {w_hi:.3f}.\n"
            "  Rebuild mosaiq_simulator and rerun the Figure 5 export."
        )

    background = load_s_ee_background(structure_factor_path, display_bounds)

    disp_rs, disp_t = parse_run_metadata(dispersion_path)
    params_path = dispersion_path
    gamma_e, tau_e = parse_electron_reduced_params(params_path)
    bare_spe = bare_imaginary_lindhard_along_roots(q, omega_p, tau_e, gamma_e)

    fig = plot_dispersion(
        q,
        omega_p,
        landau,
        bohm_gross,
        background,
        bare_spe,
        display_bounds,
        temperature_k=disp_t if disp_t is not None else MANUSCRIPT_FIG5_T_KELVIN,
    )

    saved_path = save_figure(fig, OUTPUT_FIG.stem)
    plt.close(fig)

    n_roots = int(np.sum(np.isfinite(omega_p)))
    print(f"Figure 5 data: r_s={disp_rs}, T={disp_t} K")
    print(f"  dispersion: {dispersion_path.name}")
    print(f"  structure factor: {structure_factor_path.name} (omega_max={omega_max:.3f})")
    print(f"Loaded {q.size} q points ({n_roots} finite plasmon roots).")
    if background is None:
        print(f"Warning: no background map from {structure_factor_path}; dispersion-only plot saved.")
    print(f"Saved {saved_path}")


if __name__ == "__main__":
    main()
