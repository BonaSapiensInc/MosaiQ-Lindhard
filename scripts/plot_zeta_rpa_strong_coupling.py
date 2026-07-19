#!/usr/bin/env python3
"""Plot Phase Z1 strong-coupling sweep: RPA failure counts and Experimental rescues.

Reads:
  output/zeta_rpa_strong_coupling_sweep.dat
  output/zeta_rpa_strong_coupling_rescue.dat (optional)

Writes:
  output/zeta_rpa_strong_coupling.pdf
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SWEEP = ROOT / "output" / "zeta_rpa_strong_coupling_sweep.dat"
RESCUE = ROOT / "output" / "zeta_rpa_strong_coupling_rescue.dat"
OUT = ROOT / "output" / "zeta_rpa_strong_coupling.pdf"


def main() -> int:
    if not SWEEP.is_file():
        print(f"Missing {SWEEP}; run verify_zeta_rpa_scalar first.", file=sys.stderr)
        return 1

    sweep = np.loadtxt(SWEEP, comments="#")
    if sweep.ndim == 1:
        sweep = sweep.reshape(1, -1)

    # columns: rs Gamma v_scale n_points n_rpa_fail n_prod_fail n_exp_fail n_rescue max_abs_rpa max_abs_exp
    rs = sweep[:, 0]
    gamma = sweep[:, 1]
    v_scale = sweep[:, 2]
    n_rpa_fail = sweep[:, 4]
    n_rescue = sweep[:, 7]
    max_abs_rpa = sweep[:, 8]
    max_abs_exp = sweep[:, 9]

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), constrained_layout=True)

    # Panel A: failure / rescue counts vs Gamma at fixed rs=1, aggregated over v_scale
    mask = np.isclose(rs, 1.0)
    if np.any(mask):
        g_unique = np.unique(gamma[mask])
        fail_sum = [n_rpa_fail[mask & np.isclose(gamma, g)].sum() for g in g_unique]
        rescue_sum = [n_rescue[mask & np.isclose(gamma, g)].sum() for g in g_unique]
        axes[0].plot(g_unique, fail_sum, "k-o", ms=5, label="RPA fail/extreme")
        axes[0].plot(g_unique, rescue_sum, "C3-s", ms=5, label="Experimental rescue")
        axes[0].set_xlabel(r"$\Gamma$")
        axes[0].set_ylabel("grid-point count (summed $v$)")
        axes[0].set_title(r"$r_s=1$ (incl. pole-parking $v$)")
        axes[0].legend(frameon=False, fontsize=8)
    else:
        axes[0].text(0.5, 0.5, "no rs=1 slice", ha="center", va="center")
        axes[0].set_axis_off()

    # Panel B: max |chi| RPA vs Experimental across rescues or sweep
    if RESCUE.is_file() and RESCUE.stat().st_size > 0:
        try:
            rescue = np.loadtxt(RESCUE, comments="#")
        except ValueError:
            rescue = np.empty((0, 15))
        if rescue.ndim == 1 and rescue.size > 0:
            rescue = rescue.reshape(1, -1)
        if rescue.size > 0:
            # AbsChiRPA, AbsChiExp
            sc = axes[1].scatter(
                rescue[:, 8], rescue[:, 12], c=rescue[:, 1], cmap="viridis", s=18
            )
            axes[1].plot([1e-2, 1e8], [1e-2, 1e8], "k--", lw=0.8, alpha=0.5)
            axes[1].set_xscale("log")
            axes[1].set_yscale("log")
            axes[1].set_xlabel(r"$|\chi^{\mathrm{RPA}}|$")
            axes[1].set_ylabel(r"$|\chi^{\zeta}_{\mathrm{exp}}|$")
            axes[1].set_title("Rescue points (color = Γ)")
            cb = fig.colorbar(sc, ax=axes[1], fraction=0.046)
            cb.set_label(r"$\Gamma$")
        else:
            axes[1].text(0.5, 0.5, "no rescue points", ha="center", va="center")
            axes[1].set_axis_off()
    else:
        # Fallback: max abs from sweep
        axes[1].scatter(max_abs_rpa, max_abs_exp, c=gamma, cmap="viridis", s=18)
        axes[1].set_xscale("log")
        axes[1].set_yscale("log")
        axes[1].set_xlabel(r"max $|\chi^{\mathrm{RPA}}|$")
        axes[1].set_ylabel(r"max $|\chi^{\zeta}_{\mathrm{exp}}|$")
        axes[1].set_title("Sweep maxima")

    fig.suptitle("Phase Z1: strong-coupling RPA stress vs ZetaRPA_Experimental", fontsize=11)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT)
    print(f"Wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
