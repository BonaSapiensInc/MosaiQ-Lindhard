#!/usr/bin/env bash
# Regenerate manuscript Figures 3 and 4 (run from repository root).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

SIM="./simulator/build/mosaiq_simulator"
if [[ ! -x "$SIM" ]]; then
  echo "Error: build the simulator first:" >&2
  echo "  cmake -S simulator -B simulator/build -DCMAKE_BUILD_TYPE=Release" >&2
  echo "  cmake --build simulator/build --target mosaiq_simulator" >&2
  exit 1
fi

FIG3_RS="${1:-1}"
FIG3_T="${2:-10000}"
FIG4_RS="${3:-2}"
FIG4_T="${4:-10000}"

echo "=== Figure 3: dynamic S(q,omega) at r_s=${FIG3_RS}, T=${FIG3_T} K (ee/ei q<=4; S_ii q<=50) ==="
"$SIM" "$FIG3_RS" "$FIG3_T"
cp output/output_structure_factor.dat "output/output_structure_factor_rs${FIG3_RS}.dat"
python3 scripts/plot_Sqw.py

echo "=== Figure 4: plasmon dispersion at r_s=${FIG4_RS}, T=${FIG4_T} K (q in [0.1, 50]) ==="
"$SIM" "$FIG4_RS" "$FIG4_T"
cp output/output_structure_factor.dat "output/output_structure_factor_rs${FIG4_RS}.dat"
cp output/output_plasmon_dispersion.dat "output/output_plasmon_dispersion_rs${FIG4_RS}.dat"
python3 scripts/plot_plasmon_dispersion.py

echo "Done. Figures written to output/S_*_contour.pdf and output/plasmon_dispersion.pdf"
