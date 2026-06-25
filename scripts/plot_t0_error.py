#!/usr/bin/env python3
"""Publication-quality T=0 Lindhard limit validation figure from MosaiQ-Lindhard output (PRE style)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_FILE = PROJECT_ROOT / "output" / "output_t0_limit.dat"
OUTPUT_PDF = PROJECT_ROOT / "output" / "t0_limit_validation.pdf"

COLOR_RE = "#1f3b73"
COLOR_IM = "#8b1e1e"
COLOR_ERR_RE = "#1f3b73"
COLOR_ERR_IM = "#8b1e1e"


def load_t0_limit(path: Path) -> dict[str, np.ndarray]:
    if not path.is_file():
        raise SystemExit(f"Error: data file not found: {path}")

    table = np.loadtxt(path, comments="#")
    if table.ndim != 2 or table.shape[1] < 8:
        raise SystemExit(f"Error: expected 8 columns in {path}, got shape {table.shape}")

    return {
        "q": table[:, 0],
        "omega": table[:, 1],
        "re_exact": table[:, 2],
        "im_exact": table[:, 3],
        "re_engine": table[:, 4],
        "im_engine": table[:, 5],
        "error_re": table[:, 6],
        "error_im": table[:, 7],
    }


def parse_metadata(path: Path) -> dict[str, str]:
    metadata: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            if "q_bar =" in line:
                for token in line.replace("#", "").split():
                    if "=" in token:
                        key, value = token.split("=", 1)
                        metadata[key.strip()] = value.strip()
    return metadata


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "legend.fontsize": 9,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.linewidth": 0.8,
            "grid.alpha": 0.25,
            "grid.linewidth": 0.5,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def plot_t0_validation(data: dict[str, np.ndarray], q_label: str) -> Figure:
    configure_matplotlib()

    omega = data["omega"]

    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=(6.5, 5.5),
        gridspec_kw={"height_ratios": [1.35, 1.0], "hspace": 0.08},
    )

    ax_top.plot(
        omega,
        data["re_exact"],
        color=COLOR_RE,
        linewidth=1.8,
        label=r"$\Re\tilde{\chi}_{T=0}$ (exact)",
    )
    ax_top.plot(
        omega,
        data["im_exact"],
        color=COLOR_IM,
        linewidth=1.8,
        label=r"$\Im\tilde{\chi}_{T=0}$ (exact)",
    )
    ax_top.scatter(
        omega,
        data["re_engine"],
        s=14,
        facecolors="none",
        edgecolors=COLOR_RE,
        linewidths=0.9,
        marker="o",
        label=r"Engine $\Re\tilde{\chi}$",
        zorder=3,
    )
    ax_top.scatter(
        omega,
        data["im_engine"],
        s=14,
        facecolors="none",
        edgecolors=COLOR_IM,
        linewidths=0.9,
        marker="s",
        label=r"Engine $\Im\tilde{\chi}$",
        zorder=3,
    )

    ax_top.set_ylabel(r"$\tilde{\chi}(\bar{q}, \bar{\omega})$ [$D(\epsilon_\mathrm{F})$ units]")
    ax_top.set_title(rf"Highly degenerate limit at $\bar{{q}} = {q_label}$", fontsize=11)
    ax_top.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.35)
    ax_top.legend(loc="upper right", frameon=True, framealpha=0.92, edgecolor="0.75")

    positive_re = np.maximum(data["error_re"], np.finfo(float).tiny)
    positive_im = np.maximum(data["error_im"], np.finfo(float).tiny)

    ax_bottom.plot(
        omega,
        positive_re,
        color=COLOR_ERR_RE,
        linewidth=1.6,
        label=r"$|\Delta\Re|$",
    )
    ax_bottom.plot(
        omega,
        positive_im,
        color=COLOR_ERR_IM,
        linewidth=1.6,
        linestyle="--",
        dashes=(5, 3),
        label=r"$|\Delta\Im|$",
    )

    ax_bottom.set_yscale("log")
    ax_bottom.set_ylabel(r"Absolute error")
    ax_bottom.set_xlabel(r"Reduced frequency $\bar{\omega} = \hbar\omega / \epsilon_\mathrm{F}$")
    ax_bottom.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.35)
    ax_bottom.legend(loc="upper right", frameon=True, framealpha=0.92, edgecolor="0.75")

    fig.align_ylabels([ax_top, ax_bottom])
    fig.subplots_adjust(hspace=0.08)
    return fig


def main() -> None:
    metadata = parse_metadata(DATA_FILE)
    data = load_t0_limit(DATA_FILE)

    q_values = np.unique(data["q"])
    if q_values.size != 1:
        raise SystemExit(
            f"Error: expected a single q_bar in {DATA_FILE}, found {q_values.size} values."
        )

    q_label = metadata.get("q_bar", f"{q_values[0]:g}")
    fig = plot_t0_validation(data, q_label)

    finite_re_errors = data["error_re"][np.isfinite(data["error_re"])]
    finite_im_errors = data["error_im"][np.isfinite(data["error_im"])]

    OUTPUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PDF)
    plt.close(fig)

    print(f"Loaded {data['omega'].size} omega points at q_bar = {q_values[0]:g}.")
    if finite_re_errors.size:
        print(f"max |Delta Re| = {np.max(finite_re_errors):.3e}")
    else:
        print("max |Delta Re| = n/a (no finite reference points)")
    print(f"max |Delta Im| = {np.max(finite_im_errors):.3e}")
    print(f"Saved {OUTPUT_PDF}")


if __name__ == "__main__":
    main()
